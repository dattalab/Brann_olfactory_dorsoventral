import numpy as np
import pandas as pd
from scipy import stats
from sklearn.neighbors import NearestNeighbors

from dv_score import io


def filter_bintu(df_pred_all):
    is_bad1 = (df_pred_all.rep_num == 1) & (df_pred_all.mouse_num.isin((1, 4)))
    is_bad2 = (df_pred_all.rep_num == 2) & (df_pred_all.mouse_num.isin((1, 4, 6)))
    is_bad_either = (is_bad1) | (is_bad2)
    df_pred_keep = df_pred_all[~is_bad_either].copy()
    return df_pred_keep


def rotate_section(
    df_best, x_col="absolute_center_y", y_col="absolute_center_x", rot_deg=9
):
    x = df_best[x_col].values
    y = -1 * df_best[y_col].values
    coords = np.vstack([x, y]).T
    angle = np.deg2rad(rot_deg)
    R = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    rot_coords = coords @ R
    nn = NearestNeighbors(n_neighbors=5).fit(coords)
    distances, _ = nn.kneighbors()
    isolated = distances.mean(1) > 20
    rot_keep = rot_coords[~isolated]
    return rot_keep, isolated


def get_c(_df):
    return stats.spearmanr(_df.AP, _df.ab_dist).correlation


def get_apical_basal():
    df_res_both = io.load_merfish_AP()
    res_dict = df_res_both.groupby("level_0")
    df_mature3_mean = io.load_gep_mean(to_pandas=True)
    df_bm = io.load_olfr_info()
    olfr_dv = (df_mature3_mean["DV"] >= 60).map({False: "ventral", True: "dorsal"})
    olfr_all = (
        df_bm.OR_class.loc[df_mature3_mean.index]
        + "_"
        + olfr_dv.loc[df_mature3_mean.index]
    )
    corrs = {}
    df_subsets = {}
    for k, df_res2 in res_dict:
        g_vc = df_res2.gene.value_counts()
        g_keep = g_vc.index[g_vc >= 20]
        df_res2_subset = df_res2[
            (df_res2.close_dist) & (df_res2.gene.isin(g_keep))
        ].merge(olfr_all.rename("olfr_type"), left_on="gene", right_index=True)
        df_res2_subset["is_dorsal"] = df_res2_subset.DV >= 60
        ap_mean_dist = df_res2_subset.groupby(["gene", "is_dorsal"], as_index=False)[
            ["AP", "ab_dist"]
        ].mean()

        subsets = {
            "all": np.ones_like(ap_mean_dist.is_dorsal),
            "ventral": ~ap_mean_dist.is_dorsal,
            "dorsal": ap_mean_dist.is_dorsal,
        }

        df_subsets[k] = df_res2_subset

        corrs[k] = {kk: get_c(ap_mean_dist[v]) for kk, v in subsets.items()}
    df_corr = pd.DataFrame(corrs).T
    return df_subsets, df_corr
