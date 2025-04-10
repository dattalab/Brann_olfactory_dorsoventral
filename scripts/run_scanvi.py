import argparse
import logging
import sys
from pathlib import Path

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


def run(fn, vae_folder, n_epochs):
    logging.info(fn)
    adata_fn = "Nextflow_v3_all_samples_vst_hvg_counts.h5ad"
    out_folder = fn.parents[0]
    logging.info("Loading adata_fn")
    adata_all = sc.read(adata_fn)
    logging.info("loaded.")
    to_query = read_list(fn)
    logging.info("Running for sources:")
    logging.info(to_query)
    is_in_query = adata_all.obs.source.isin(to_query)
    adata_query = adata_all[is_in_query, :].copy()
    del adata_all
    logging.info(adata_query)
    adata_query.obs["celltype"] = "Unknown"
    vae_q = scvi.model.SCANVI.load_query_data(
        adata_query,
        str(vae_folder),
    )
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
    logging.info("Saved")


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
    run(args.fn, args.vae, args.epochs)


if __name__ == "__main__":
    main()
