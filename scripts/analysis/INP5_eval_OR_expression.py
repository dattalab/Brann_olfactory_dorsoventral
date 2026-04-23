import numpy as np
import pandas as pd
import polars as pl

from dv_score import io

adata = ""
save_data_folder = ""

df_obs = io.load_development("cell_obs")
df_olfr = io.load_development("olfr")
df_olfr_mat = io.load_development("olfr_mat")


# wide to long dataframe
X_or = adata[:, df_olfr_mat.columns].X.tocoo()
has_olfr = X_or.data >= 1
df_pl3 = pl.DataFrame(
    dict(row=X_or.row[has_olfr], col=X_or.col[has_olfr], expr=X_or.data[has_olfr])
)

# info for each cell
df_row_meta = pl.DataFrame(
    df_obs[["pseudotime_rank", "n_OR_1", "n_OR_3"]]
    .join(df_olfr.olfr_max)
    .reset_index(drop=True)
    .rename_axis(index="row")
    .reset_index()
)

# mapping for each OR
df_col_meta = pl.DataFrame(
    data=[np.arange(len(df_olfr_mat.columns)), list(df_olfr_mat.columns)],
    schema=["col", "top_Olfr"],
)

# DV scores for each OR
df_mature = io.load_df_dv(apply_filter=True)
df_mill_dv_mean = df_mature.group_by("top_Olfr").agg(pl.col("DV").mean())
df_mill_dv_mean_pd = df_mill_dv_mean.to_pandas().set_index("top_Olfr")
bm_pl = pl.DataFrame(
    io.load_olfr_info()["OR_class"].rename_axis("top_Olfr").reset_index()
)
df_or_mean = (
    df_mill_dv_mean.with_columns(pl.col("top_Olfr").cast(str))
    .join(bm_pl, on="top_Olfr")
    .join(df_col_meta, on="top_Olfr")
    .with_columns(
        pl.col("DV").qcut(10, labels=(np.arange(10) + 1).astype(str)).alias("cat")
    )
)

# join with pseudotime of each cell and DV score of each OR
df_pl3_pseudotime = (
    df_pl3.join(df_row_meta, on="row")
    .join(df_or_mean, on="col", how="inner")
    .with_columns(has_OR=pl.col("n_OR_3") >= 1, above_thresh=pl.col("expr") >= 3)
)

# cells to eval
df_pl3_c2 = df_pl3_pseudotime.filter(
    (pl.col("pseudotime_rank") > 0.4) & (pl.col("OR_class") == "class2")
)
df_pl3_keep = df_pl3_c2.filter((pl.col("has_OR")) & (pl.col("above_thresh")))
df_pl3_low = df_pl3_c2.filter(pl.col("expr") < 3)


n_first = 30
olfr_keep = (
    df_pl3_keep["top_Olfr"]
    .value_counts()
    .filter(pl.col("count") >= n_first)["top_Olfr"]
    .unique()
)

cat2 = pd.qcut(df_mill_dv_mean_pd.loc[olfr_keep].DV, 10)
new_codes = pl.DataFrame((cat2.cat.codes + 1).rename("cat2").reset_index())
print(len(olfr_keep))
df_pl_sort = (
    df_pl3_keep.filter(pl.col("top_Olfr").is_in(olfr_keep))
    .join(new_codes, on="top_Olfr")
    .sort(["top_Olfr", "pseudotime_rank"])
)

# earliest cells expressing each OR
df_first_mean = df_pl_sort.group_by(["cat2", "top_Olfr"]).agg(
    pl.col("pseudotime_rank").quantile(0.1)
)

# quantiles of pre-choice cells
_df = df_pl3_low.to_pandas()
breaks = [l.left for l in cat2.cat.categories]
breaks[0] = -np.inf
breaks.append(np.inf)
c = pd.cut(_df["DV"], breaks, include_lowest=True)
_df["cat2"] = c.cat.codes + 1
gb = _df.groupby(["cat2", "top_Olfr"], observed=True)
quants = {
    f"pre-{q}": gb.pseudotime_rank.quantile(q).reset_index() for q in [0.25, 0.5, 0.75]
}
quants["post-choice"] = df_first_mean.to_pandas()
to_plot = pd.concat(quants).reset_index()


df_export = df_pl3_c2.to_pandas().copy()
df_export2 = df_export.rename(
    columns={
        "row": "cell_num",
        "col": "olfr_num",
        "expr": "olfr_expression",
        "n_OR_1": "n_OR_1_umi",
        "n_OR_3": "n_OR_3_umi",
        "olfr_max": "cell_olfr_max",
        "has_OR": "cell_any_OR_3_umi",
    }
).drop(columns=["cat", "OR_class"])
df_export2.to_parquet(save_data_folder / "INP_cells_OR_expression_long.parquet")
df_pl3_pseudotime.write_parquet(
    save_data_folder / "INP5_OR_expression_explode_pre_post_choice.parquet"
)
to_plot.to_parquet(
    save_data_folder / "INP5_pseudotime_rank_vs_dv_score_bin_pre_post_choice.parquet"
)


###############
# weight OR expression by DV score
df_gep_mean = io.load_gep_mean().to_pandas().set_index('top_Olfr')
olfr_names = list(df_gep_mean.index)
df_olfr_mat = io.load_development("olfr_mat")

olfr_mat = df_olfr_mat.select(pl.col(olfr_names)).to_numpy()
n_cells, n_olfr = olfr_mat.shape

olfr_max = olfr_mat.max(1)
olfr_argsort = np.argsort(olfr_mat, axis=1)
olfr_sort = np.sort(olfr_mat, axis=1)
has_any = (olfr_sort > 0).any(0)
# DV score for each OR

olfr_idx = np.empty(n_olfr)
olfr_idx.fill(np.nan)

dv_means = df_gep_mean.DV
in_dv = df_olfr_mat.columns.isin(dv_means.index)
olfr_idx[in_dv] = dv_means.loc[olfr_names[in_dv]]
dv_mat = olfr_idx[olfr_argsort[:, has_any]]
dv_mat[olfr_sort[:, has_any] == 0] = np.nan

# weighted-mean, weight DV score by expression
weights = olfr_sort[:, has_any] / olfr_sort.sum(1)[:, None]
weighted_dv = np.nansum(dv_mat * weights, 1)
