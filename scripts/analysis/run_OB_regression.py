from warnings import simplefilter

import numpy as np
import pandas as pd
import sklearn.linear_model as lm
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from tqdm.auto import trange

from dv_score import io

simplefilter("ignore", category=ConvergenceWarning)

df_hvg_new = io.load_hvg_mean()
df_geom_dv = io.load_glomeruli()
df_median = df_geom_dv.groupby(["domain", "top_Olfr"], as_index=False).median(
    numeric_only=True
)

# constants
GEP_dict = {
    "DV": ["Dorsal", "Ventral", "DV"],
    "AP": ["Anterior", "Posterior", "AP"],
    "ES": ["High", "Low", "ES"],
}
GEP_dict["DV_AP"] = GEP_dict["DV"] + GEP_dict["AP"]
pos_cols = ["x_flipped", "y", "z"]
n_rep = 1000

lr_en = lm.ElasticNet(alpha=1, l1_ratio=0.9)
lr_steps = [("scale", StandardScaler()), ("mdl", lr_en)]
lr_pipe = Pipeline(lr_steps)

all_preds = {}
for d, df_med2 in df_median.groupby("domain"):
    print(d)
    # can fit each Y column at once
    Y_test = df_med2[pos_cols]
    col_names = pos_cols + ["all"]
    preds = {}
    for k, this_gep in GEP_dict.items():
        X_test = df_med2[this_gep]
        out_mat = np.zeros((n_rep, 4))
        for i in trange(n_rep, ncols=80):
            kf = KFold(n_splits=5, shuffle=True, random_state=i)
            splits = list(kf.split(X_test))
            preds1 = cross_val_predict(lr_pipe, X_test, Y_test, cv=splits)
            mse = (Y_test - preds1) ** 2
            out_mat[i, :3] = np.median(np.sqrt(mse), axis=0)
            out_mat[i, 3] = np.median(np.sqrt(mse.sum(1)), axis=0)
        preds[k] = pd.DataFrame(out_mat, columns=col_names)
    all_preds[d] = pd.concat(preds)

# regression with variable genes
all_preds_hvg = {}
for d, df_med2 in df_median.groupby("domain"):
    print(d)
    Y_test = df_med2[pos_cols]
    col_names = pos_cols + ["all"]
    preds = {}
    olfrs = df_med2.index.get_level_values(1)
    X_gene = df_hvg_new.loc[olfrs]
    out_mat = np.zeros((n_rep, 4))
    for i in trange(n_rep, ncols=80):
        kf = KFold(n_splits=5, shuffle=True, random_state=i)
        splits = list(kf.split(X_gene))
        preds1 = np.zeros_like(Y_test)
        for j, col in enumerate(Y_test.columns):
            preds1[:, j] = cross_val_predict(lr_pipe, X_gene, Y_test[col], cv=splits)
        mse = (Y_test - preds1) ** 2
        out_mat[i, :3] = np.median(np.sqrt(mse), axis=0)
        out_mat[i, 3] = np.median(np.sqrt(mse.sum(1)), axis=0)
        preds["HVG"] = pd.DataFrame(out_mat, columns=col_names)
    all_preds_hvg[d] = pd.concat(preds)

df_mae_predictions = pd.concat(
    {"GEP": pd.concat(all_preds), "HVG": pd.concat(all_preds_hvg)}
)
