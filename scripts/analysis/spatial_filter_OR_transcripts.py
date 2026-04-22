import numpy as np
import pandas as pd
from sklearn import metrics
from tqdm import tqdm

from dv_score import io

THRESH = 10
N_LOCAL = 3

# codebook and detected transcripts on GEO at GSE297068
df_codebook = pd.read_csv("codebook.csv")
# for an example region
REGION_NAME = "region_3"
df_tr = pd.read_csv(f"{REGION_NAME}_detected_transcripts.csv")
# rotate
theta = np.deg2rad(30.2 + 90)
ct = np.cos(theta)
st = np.sin(theta)
mat = np.array([[ct, -1 * st], [st, ct]])
rot_coord = df_tr[["global_x", "global_y"]].values @ mat
df_tr[["gx", "gy"]] = rot_coord

df_olfr_codes = df_codebook[
    (df_codebook.barcodeType == "merfish")
    & (df_codebook.name.str.contains("Olfr|Taar"))
]
olfr_set = set(df_olfr_codes.name)
is_olfr = df_tr.gene.isin(olfr_set)
df_olfr = df_tr[is_olfr].copy()

df_mill_dv_mean = io.load_gep_mean(to_pandas=True)
df_olfr_merge = df_olfr.merge(df_mill_dv_mean, left_on="gene", right_index=True)

# only keep "cells" with enough OR transcripts within local distance
ind_keep = pd.Index([])
outs = []
for _, _df in tqdm(df_olfr_merge.groupby("gene")):
    dmat = metrics.pairwise_distances(_df)
    local_sum = (dmat <= THRESH).sum(1)
    has_local = local_sum >= N_LOCAL + 1
    ind_keep = ind_keep.append(_df.index[has_local])
    outs.append(pd.DataFrame(local_sum, index=_df.index, columns=["n_neigh"]))

df_olfr_neigh = df_olfr_merge.join(pd.concat(outs))
df_olfr_local = df_olfr_neigh.loc[ind_keep]
df_olfr_local.to_parquet(f"{REGION_NAME}_df_olfr_local.parquet")
