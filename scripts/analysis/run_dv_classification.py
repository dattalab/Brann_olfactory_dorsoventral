import argparse
import itertools
import time
import uuid
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.svm import SVC

from dv_score.util import get_cores


def make_pipe():
    pipe = Pipeline(
        [
            ("scale", StandardScaler()),
            ("clf", SVC(C=0.5, kernel="rbf")),
        ]
    )
    return pipe


def parse_args():
    path_type = lambda p: Path(p).absolute()
    # parse args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-t",
        dest="thresh",
        help="number of cells per OR threshold",
        type=int,
        default=150,
    )
    parser.add_argument(
        "-r",
        dest="repeats",
        help="number of repeats of pipeline",
        type=int,
        default=10,
    )
    parser.add_argument(
        "-p", dest="pick", help="number of cells per OR", type=int, default=150
    )
    parser.add_argument("-n", dest="run_index", help="restart", type=int, default=25)
    parser.add_argument(
        "--data",
        dest="data",
        type=path_type,
        default="DV_score_all_mature_OSNs.parquet",
    )
    parser.add_argument(
        "--out",
        dest="out_folder",
        type=path_type,
        default="/mnt/disks/sec/David/MOE/development/data/mature_dv_classification",
    )
    return parser.parse_args()


def run(DATA, THRESH, N_SUBSET, RUN_INDEX, OUT_FOLDER):
    """
    Run the DV classification pipeline.

    Args:
        DATA (str): Path to the data file.
        THRESH (int): Threshold value for number of cells per OR.
        N_SUBSET (int): Number of cells per OR to subsample.
        RUN_INDEX (int): Run index.
        OUT_FOLDER (str): Path to the output folder.

    Returns:
        None
    """
    df_mature3 = pd.read_parquet(
        DATA, columns=["Dorsal", "Ventral", "DV", "top_Olfr", "has_OR"]
    )
    df_mature3_use = df_mature3[
        (df_mature3.has_OR) & (~df_mature3.index.str.contains("opto"))
    ]
    vc_all = df_mature3_use.top_Olfr.value_counts()
    at_least_n = vc_all.index[vc_all >= THRESH]
    df_mature2_keep = df_mature3_use[df_mature3_use.top_Olfr.isin(at_least_n)]
    le = LabelEncoder().fit(df_mature2_keep.top_Olfr)
    idx = le.transform(df_mature2_keep.top_Olfr)
    dv_score = df_mature2_keep.DV.values
    uq_idx = np.unique(idx)
    N_OR = len(uq_idx)
    np.random.seed(RUN_INDEX)
    rands = []
    for i in uq_idx:
        is_i = np.where(idx == i)[0]
        rands.append(np.random.choice(is_i, (N_SUBSET), replace=True))
    rands = np.hstack(rands)

    y_subset = idx[rands]
    X_subset = df_mature2_keep[["Dorsal", "Ventral", "DV"]].values[rands]
    dv_rand = dv_score[rands]
    print(dv_rand.shape)
    uq_y = np.unique(y_subset)
    combos = list(itertools.combinations(uq_y, 2))
    PIPE = make_pipe()
    print(
        f"Run {RUN_INDEX}: {len(combos)} combos with thresh {THRESH} and {N_SUBSET} samples per OR."
    )

    def run_combo(combo, N_SPLIT=5):
        is_combo = np.in1d(y_subset, combo)
        y_combo = y_subset[is_combo]
        X_combo = X_subset[is_combo]
        skf = StratifiedKFold(n_splits=N_SPLIT, shuffle=True, random_state=RUN_INDEX)
        pred = cross_val_predict(PIPE, X_combo, y_combo, cv=skf, n_jobs=1)
        obs = (pred == y_combo).mean()
        y_shuff = np.random.permutation(y_combo)
        pred = cross_val_predict(PIPE, X_combo, y_combo, cv=skf, n_jobs=1)
        shuff = (pred == y_shuff).mean()

        is_a = y_combo == combo[0]
        X_dv_combo = dv_rand[is_combo]
        dv_perm = np.random.permutation(X_dv_combo)
        frac_overlap = ((X_dv_combo[is_a, None] > X_dv_combo[~is_a])).mean()
        shuff_overlap = ((dv_perm[is_a, None] > dv_perm[~is_a])).mean()
        is_b = ~is_a if X_dv_combo[is_a].mean() < X_dv_combo[~is_a].mean() else is_a
        auroc = metrics.roc_auc_score(is_b, X_dv_combo)
        # afterwards correct if less than 0.5
        auroc_shuff = metrics.roc_auc_score(is_b, dv_perm)
        return (obs, shuff, frac_overlap, shuff_overlap, auroc, auroc_shuff)

    cores = (get_cores() - 1) // 5 * 5
    res = Parallel(n_jobs=cores, verbose=True)(
        delayed(run_combo)(combo) for combo in combos
    )
    df_res = pd.DataFrame(
        res,
        columns=[
            "obs",
            "shuff",
            "frac_overlap",
            "frac_overlap_shuff",
            "auroc",
            "auroc_shuff",
        ],
        index=combos,
    )
    this_uuid = uuid.uuid4().hex
    to_save = OUT_FOLDER / f"{THRESH}_{N_SUBSET}_{N_OR}"
    to_save.mkdir(exist_ok=True, parents=True)
    outfn = (
        to_save
        / f"DV_pairwise_classification_{THRESH}_{N_SUBSET}_{N_OR}_{RUN_INDEX}_{this_uuid}.parquet"
    )
    df_res.to_parquet(outfn)
    print(f"Saved to {outfn}")


def main(args):
    N_REP = args.repeats
    start_time = time.time()
    for i in range(N_REP):
        print(f"Starting run {i + 1} of {N_REP}")
        run(args.data, args.thresh, args.pick, args.run_index, args.out_folder)
        elapsed = (time.time() - start_time) / 60
        print(f"Finished run {i + 1} after {elapsed:.4f} minutes.")


if __name__ == "__main__":
    args = parse_args()
    print(args)
    main(args)
