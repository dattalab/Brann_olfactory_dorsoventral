import numpy as np
import numpy_groupies as npg
import pandas as pd
from scipy.spatial import distance
from sklearn import cluster
from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import LabelEncoder
from tqdm import tqdm

from dv_score import io

# global aligned coordinates
pos_cols = ["gx", "gy"]

# see transcript files on GEO GSE297068
df_tr_all = pd.read_parquet("region_0_3_combined_detected_transcripts.parquet")
# KNN-based clustering using non-OR transcripts ala Baysor
df_labs = pd.read_parquet("Coronal2_region_0_3_knn_labels.parquet")

regions = df_tr_all.level_0.unique()
df_locs = {}
for r in regions:
    df_locs[r] = pd.read_parquet(f"{r}_df_olfr_local.parquet")
df_locs_all = pd.concat(df_locs).reset_index()

df_locs_merge = df_locs_all.merge(
    df_tr_all[["level_0", "level_1", "gx", "gy"]],
    on=["level_0", "level_1"],
    suffixes=("_r", None),
)

df_mature3_mean = io.load_gep_mean(to_pandas=True)
df_locs_merge = df_locs_merge.merge(
    df_mature3_mean[["DV", "AP", "ES"]],
    how="left",
    left_on="gene",
    right_index=True,
    suffixes=("_old", None),
)

dfs = {}
db = cluster.DBSCAN(eps=5, min_samples=3)
for olfr, _df in tqdm(df_locs_merge.groupby("gene")):
    db.fit(_df[["global_x", "global_y"]])
    _df["cell_labels"] = db.labels_
    _df["olfr_cell_labels"] = [olfr + "_" + l for l in db.labels_.astype(str)]
    dfs[olfr] = _df

df_cat = pd.concat(dfs)
df_cat.index.names = ["olfr", "level_2"]
df_olfr_labels = df_cat.reset_index()
df_olfr_cells = df_olfr_labels[df_olfr_labels.cell_labels != -1]


df_olfr_cells_labs = df_olfr_cells.merge(
    df_labs, left_on="level_1", right_index=True, how="left"
)
cell_clust_counts = pd.crosstab(
    df_olfr_cells_labs.olfr_cell_labels, df_olfr_cells_labs.osn_labs
)

min_tr_per_cell = df_olfr_cells.olfr_cell_labels.value_counts().min()

# keep only cells with enough OSN transcripts
good_cells = cell_clust_counts.index[cell_clust_counts[0] >= min_tr_per_cell]

df_olfr_cells = df_olfr_cells_labs[
    (df_olfr_cells_labs.olfr_cell_labels.isin(good_cells))
    & (df_olfr_cells_labs.osn_labs == 0)
]

# look at pairwise distances between transcripts within cell
nc = df_olfr_cells.olfr_cell_labels.nunique()
maxes = np.zeros(nc)
meds = np.zeros(nc)
for i, (_, _df) in tqdm(enumerate(df_olfr_cells.groupby("olfr_cell_labels"))):
    dm = distance.pdist(_df[pos_cols])
    maxes[i] = dm.max()
    meds[i] = np.median(dm)
bad_cells = (meds >= 10) | (maxes >= 30)


BAD_ORS = ["Olfr672", "Olfr1225", "Olfr640", "Olfr1065"]  # in autofluoresence patches
CELL_THRESH = 5

uq_cell_names = np.unique(df_olfr_cells.olfr_cell_labels)
cells_keep = uq_cell_names[~bad_cells]
is_good_olfr = ~df_olfr_cells.gene.isin(BAD_ORS)
is_good_cell = df_olfr_cells.olfr_cell_labels.isin(cells_keep)

n_cells_per_olfr = df_olfr_cells.groupby("gene").olfr_cell_labels.nunique()
is_good_or = df_olfr_cells.gene.isin(
    n_cells_per_olfr.index[n_cells_per_olfr >= CELL_THRESH]
)

df_olfr_cells = df_olfr_cells[(is_good_olfr) & (is_good_cell) & (is_good_or)]


# KNN regression
df_cell_mean = df_olfr_cells.groupby(
    ["gene", "olfr_cell_labels"], as_index=False
).mean()
cell_mean_mean = df_cell_mean.groupby("gene").mean()
ys = df_cell_mean.DV.values
X_cells = df_cell_mean[pos_cols].values
olfr_le = LabelEncoder().fit(df_cell_mean.gene)
olfr_labels = olfr_le.transform(df_cell_mean.gene)

knn = KNeighborsRegressor(n_neighbors=10, weights="distance")
# leave out one OSN subtype at a time
loo = LeaveOneGroupOut()
preds = cross_val_predict(
    knn, X_cells, ys, cv=loo, n_jobs=-1, verbose=True, groups=olfr_labels
)
pred_mean = npg.aggregate(olfr_labels, preds, func="mean")
df_preds = pd.DataFrame(
    [preds, ys],
    index=["pred_KNN", "y_DV"],
    columns=df_cell_mean.index,
).T.join(df_cell_mean[["gene", "olfr_cell_labels", "global_x", "global_y", "gx", "gy"]])
