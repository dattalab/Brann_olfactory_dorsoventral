import numpy as np
import pandas as pd
from scipy import stats
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import LabelEncoder
from tqdm import trange

from dv_score import util


def pair_corr(df_multi, n_boot=10_000):
    clone_idx = LabelEncoder().fit_transform(df_multi.mouse_BC12)
    n_c2 = df_multi.mouse_BC12.nunique()
    clone_two = util.subset_source_OR_indices(clone_idx, 2, n_boot=n_boot)
    c2_dvm = df_multi.DV_M.values
    c2_pairs = c2_dvm[clone_two]
    pair_corrs = np.zeros((n_boot, 2))
    a = c2_pairs[::2]
    b = c2_pairs[1::2]
    rand_cells = np.stack(
        [np.random.choice(len(c2_dvm), n_c2 * 2, replace=False) for _ in range(n_boot)]
    )
    rand_dvs = c2_dvm[rand_cells]
    for i in trange(n_boot, ncols=80):
        pair_corrs[i, 0] = stats.spearmanr(a[:, i], b[:, i])[0]
        pair_corrs[i, 1] = stats.spearmanr(rand_dvs[i, ::2], rand_dvs[i, 1::2])[0]
    return pd.DataFrame(pair_corrs, columns=["observed", "shuffled"])


def within_between_clone(df_multi, n_rand=1_000):
    # used n_rnad = 10_000 for figures
    # pairwise difference between cells
    dv_dist_mat = np.abs(df_multi.DV_M.values[:, None] - df_multi.DV_M.values[None, :])
    dv_pct_mat = np.abs(df_multi.pct.values[:, None] - df_multi.pct.values[None, :])
    tri_idx = np.triu_indices(dv_dist_mat.shape[0], 1)

    df_diag_dist = util.melt_diag(
        pd.DataFrame(dv_dist_mat, index=df_multi.cell, columns=df_multi.cell),
        k=1,
    )
    assert np.allclose(df_diag_dist.value, dv_dist_mat[tri_idx])
    df_diag_dist["pct"] = dv_pct_mat[tri_idx]

    df_diag_annot = df_diag_dist.merge(
        df_multi[["cell", "BC12", "source"]], left_on="cell_1", right_on="cell"
    ).merge(df_multi[["cell", "BC12", "source"]], left_on="cell_2", right_on="cell")
    same_BC = df_diag_annot.BC12_x == df_diag_annot.BC12_y
    df_diag_annot["same_BC"] = same_BC
    same_source = df_diag_annot.source_x == df_diag_annot.source_y
    within_between = (
        df_diag_annot[(same_source)]
        .groupby(["same_BC", "source_x"], observed=True)
        .pct.median()
        .reset_index()
    )
    df_same = df_diag_annot[same_BC]
    n_same = same_BC.sum()
    rands_same = np.random.choice(df_diag_annot.shape[0], (n_same, n_rand))
    big_pct_mat = df_diag_annot.pct.values[rands_same]
    xsz = np.arange(101)
    cmat_pct = np.zeros((len(xsz), n_rand))
    for i, _x in enumerate(xsz):
        cmat_pct[i] = (big_pct_mat <= _x).mean(0)

    return within_between, df_same, cmat_pct
