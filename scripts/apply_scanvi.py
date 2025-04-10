import argparse
import gc
import logging
import pickle
import sys
from pathlib import Path

import numba
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

logging.basicConfig(
    stream=sys.stdout,
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO,
)

scvi.settings.verbosity = 30


def read_list(fn):
    with open(fn, "r") as f:
        vars = [l.rstrip("\n").strip() for l in f.readlines()]
    return vars


@numba.njit
def weighted_prediction(weights, ref_cats, N_cats):
    """Get highest weight category."""
    N = len(weights)
    predictions = np.zeros((N,), dtype=ref_cats.dtype)
    uncertainty = np.zeros((N,))
    probabilities = np.zeros((N, N_cats))
    for i in range(N):
        obs_weights = weights[i]
        obs_cats = ref_cats[i]
        best_prob = 0
        for c in np.unique(obs_cats):
            cand_prob = np.sum(obs_weights[obs_cats == c])
            probabilities[i, c] = cand_prob
            if cand_prob > best_prob:
                best_prob = cand_prob
                predictions[i] = c
                uncertainty[i] = max(1 - best_prob, 0)

    return predictions, uncertainty, probabilities


def predict_cell_types(adata, k_neigh=100):
    dev_folder = Path()
    pkl_fn = (
        dev_folder
        / "NN_index_all_samples_vst_hvg_counts_UMAP_k500_no_obsp_celltype_knn_100_neighbors.pkl"
    )
    with open(pkl_fn, "rb") as f:
        ref_nn_index = pickle.load(f)
    query_emb = adata.obsm["X_scANVI"]
    ref_neighbors, ref_distances = ref_nn_index.query(query_emb, k=k_neigh)
    df_obs = pd.read_parquet(
        dev_folder
        / "Nextflow_v3_all_samples_vst_hvg_counts_UMAP_k500_no_obsp_celltype_knn_obs_df.parquet",
        columns=["celltype_pred"],
    )
    stds = np.std(ref_distances, axis=1)
    stds = (2.0 / stds) ** 2
    stds = stds.reshape(-1, 1)
    ref_distances_tilda = np.exp(-np.true_divide(ref_distances, stds))
    weights = ref_distances_tilda / np.sum(ref_distances_tilda, axis=1, keepdims=True)
    # convert distances to affinities
    celltypes = df_obs["celltype_pred"]
    N_cats = len(np.unique(celltypes))
    ref_cats = celltypes.cat.codes.to_numpy()[ref_neighbors]
    p, u, prob = weighted_prediction(weights, ref_cats, N_cats)
    p_names = np.asarray(celltypes.cat.categories)[p]
    df_prob = pd.DataFrame(prob, index=adata.obs_names, columns=np.unique(celltypes))

    adata.obs["celltype_pred"] = p_names
    adata.obs["celltype_uncertainty"] = u
    adata.obsm["celltype_prob"] = df_prob

    sc.pp.neighbors(adata, n_neighbors=25, use_rep="X_scANVI")
    sc.tl.umap(adata)
    return adata


def run(fn, vae_folder, n_epochs, next_folder, MITO_PCT=50):
    logging.info(fn)
    out_folder = fn.parents[0]
    for f in ("adatas", "fitted_models", "embedding"):
        (out_folder / f).mkdir(exist_ok=True)
    to_query = read_list(fn)
    logging.info("Running for sources:")
    logging.info(to_query)

    prefixes = read_list(fn)
    adata_fns = [
        a
        for a in next_folder.glob("*/*.h5ad")
        if a.stem.split("_counts")[0] in prefixes
    ]
    assert set([a.parent.name for a in adata_fns]) == set(prefixes)

    adatas = []
    for f in adata_fns:
        this_adata = sc.read(f)
        adatas.append(this_adata[this_adata.obs.pct_counts_mito <= MITO_PCT].copy())

    adata_query = adatas[0].concatenate(adatas[1:], index_unique=None, join="outer")
    del adatas
    del adata_query.var

    gc.collect()
    logging.info(adata_query)
    adata_query.obs["source"] = adata_query.obs["orig_ident"].copy()
    adata_query.obs["celltype"] = "Unknown"

    scvi.model.SCANVI.prepare_query_anndata(adata_query, str(vae_folder))

    vae_q = scvi.model.SCANVI.load_query_data(
        adata_query,
        str(vae_folder),
    )
    # no longer necessary?
    vae_q._unlabeled_indices = np.arange(adata_query.n_obs)
    vae_q._labeled_indices = []
    logging.info(vae_q)
    train_kwargs_surgery = {
        "early_stopping": True,
        "early_stopping_monitor": "elbo_train",
        "early_stopping_patience": 10,
        "early_stopping_min_delta": 0.001,
        "plan_kwargs": {"weight_decay": 0.0},
    }
    vae_q.train(max_epochs=n_epochs, **train_kwargs_surgery)
    logging.info("Saving")
    name_base = fn.stem
    vae_q.save(out_folder / f"fitted_models/{name_base}", overwrite=True)

    df_embedding = pd.DataFrame(
        vae_q.get_latent_representation(), index=adata_query.obs_names
    )
    df_embedding.columns = "latent_" + df_embedding.columns.astype(str)
    df_embedding.to_parquet(out_folder / f"embedding/{name_base}.parquet")

    adata_query.obsm["X_scANVI"] = vae_q.get_latent_representation()
    adata_query.obs["predictions"] = vae_q.predict()

    # predict cell types
    logging.info("getting celltype predictions")
    adata_query = predict_cell_types(adata_query)

    adata_query.write_h5ad(out_folder / f"adatas/{name_base}.h5ad")
    logging.info("Saved anndata object")


def parse_args():
    # parse args
    parser = argparse.ArgumentParser()
    path_type = lambda p: Path(p).absolute()
    parser.add_argument("-f", "--fn", dest="fn", help="filename", type=path_type)
    parser.add_argument(
        "-v",
        "--vae",
        dest="vae",
        help="scANVI vae folder",
        type=path_type,
    )
    parser.add_argument(
        "-n",
        "--next",
        dest="next_folder",
        help="nextflow output folder",
        type=path_type,
    )
    parser.add_argument(
        "-e",
        "--epochs",
        dest="epochs",
        help="Number of epochs to train model",
        type=int,
        default=200,
    )
    return parser.parse_args()


def main():
    args = parse_args()
    run(args.fn, args.vae, args.epochs, args.next_folder)


if __name__ == "__main__":
    main()
