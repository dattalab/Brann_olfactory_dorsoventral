import numpy as np
import numpy_groupies as npg
import pandas as pd
from sklearn.preprocessing import LabelEncoder, QuantileTransformer
from tqdm import trange

from dv_score import constants, io


def delta_or(
    df_wt_dv, starts=np.arange(-400, 251, step=10), bin_width=100, only_c2=True
):
    if only_c2:
        is_c2 = (df_wt_dv.OR_class == "class2").values
    else:
        is_c2 = np.ones(len(df_wt_dv), type=bool)

    geno_idx = LabelEncoder().fit_transform(df_wt_dv.geno2)
    n_per_id = npg.aggregate(geno_idx, 1)
    dvm = df_wt_dv.DV_M.values
    in_region = (
        (dvm >= starts[:, None]) & (dvm <= starts[:, None] + bin_width) & (is_c2)
    )
    sums = npg.aggregate(geno_idx, in_region.T, axis=0, func="sum")
    sum_norm = sums / n_per_id[:, None] * 100
    df_norm_plot = pd.DataFrame(
        sum_norm, index=np.unique(df_wt_dv.geno2), columns=starts + bin_width // 2
    )
    # mouse level
    mouse_idx = LabelEncoder().fit_transform(df_wt_dv.orig_ident)
    uq_mice = np.unique(df_wt_dv.orig_ident)
    n_per_mouse = npg.aggregate(mouse_idx, 1)
    mouse_sums = npg.aggregate(mouse_idx, in_region.T, axis=0, func="sum")
    mouse_sum_norm = mouse_sums / n_per_mouse[:, None] * 100
    df_mouse_sum = pd.DataFrame(mouse_sum_norm, index=uq_mice, columns=starts)

    return df_norm_plot, df_mouse_sum


def get_norm(df, edges, cond_col):
    df["all_pct_bin"] = pd.cut(df.DV2, edges, include_lowest=True).cat.codes
    bin_ct = pd.crosstab(df.all_pct_bin, df[cond_col])
    bin_ct_norm = bin_ct / bin_ct.sum(0) * 100
    return bin_ct_norm


def log_fc(
    df_wt_dv, control_map, edges=np.arange(101, step=10), N_REP=100, cond_col="geno2"
):
    df_qt_mean = io.load_gep_mean().to_pandas().set_index("top_Olfr")[constants.GEP3]
    # renorm means from 0-1
    df_qt_mean["DV2"] = (
        QuantileTransformer().fit_transform(df_qt_mean.DV.values[:, None]) * 100
    )
    df_wt_dv_qt = df_wt_dv.merge(
        df_qt_mean.DV2,
        left_on="top_Olfr",
        right_index=True,
        how="inner",
        suffixes=(None, "_qt"),
    )
    df_wt_dv_qt["all_pct_bin"] = pd.cut(
        df_wt_dv_qt.DV2, edges, include_lowest=True
    ).cat.codes
    bin_ct = pd.crosstab(df_wt_dv_qt.all_pct_bin, df_wt_dv_qt[cond_col])
    bin_ct_norm = bin_ct / bin_ct.sum(0) * 100

    # bootstrap
    norms = []
    for i in trange(N_REP, ncols=80):
        df_sample = df_wt_dv_qt.groupby(cond_col).sample(
            frac=1, replace=True, random_state=i
        )
        norms.append(get_norm(df_sample, edges, cond_col))

    ci_dict = {}
    for g, control in control_map.items():
        norm_cols = norms[0].columns
        norm_stack = np.stack(norms)

        g_all = norm_stack[:, :, norm_cols == g].squeeze()
        c_all = norm_stack[:, :, norm_cols == control].squeeze()

        pcts = (g_all[:, None] / c_all[None]).reshape(-1, 10)
        ci_dict[g] = np.percentile(np.log2(pcts), [2.5, 97.5], axis=0)

    return bin_ct_norm, ci_dict
