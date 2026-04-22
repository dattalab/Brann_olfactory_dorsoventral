import os

import numpy as np
import numpy_groupies as npg
import pandas as pd
import polars as pl
from sklearn.preprocessing import LabelEncoder, QuantileTransformer

from dv_score.constants import GEP3, GEP3_qt


def get_cores():
    try:
        cores = int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
    except KeyError as ke:
        print("Not on SLURM, setting cores to max cpus")
        cores = os.cpu_count()
    print(f"Found {cores} cores.")
    return cores


def get_ci_vect_vectorized(x, n_boots=1000, n_samp=None, function=np.mean, pct=5):
    if isinstance(x, pd.core.series.Series):
        x = x.values
    pct /= 2
    n_vals = len(x)
    if n_samp is None:
        n_samp = n_vals
    boots = function(x[np.random.choice(n_vals, size=(n_samp, n_boots))], axis=0)
    return np.percentile(boots, [pct, 100 - pct])


def get_ci_mat_vectorized(x, n_boots=1000, n_samp=None, function=np.mean, pct=5):
    """CI for matrix of observations by variables
    * currently samples"""
    pct /= 2
    n_vals = x.shape[0]
    if n_samp is None:
        n_samp = n_vals
    boots = function(x[np.random.choice(n_vals, size=(n_samp, n_boots))], axis=0)
    return np.percentile(boots, [pct, 100 - pct], axis=0)


def pd_qt(df, geps=GEP3):
    df_qt = pd.DataFrame(
        QuantileTransformer(random_state=42).fit_transform(df[geps]) * 100,
        index=df.index,
        columns=GEP3_qt,
    )
    return df.join(df_qt, rsuffix="_qt")


def pl_qt(df: pl.DataFrame, geps=GEP3) -> pl.DataFrame:
    gep_qt = [g + "_qt" for g in geps]
    X_qt = QuantileTransformer(random_state=42).fit_transform(df.select(geps).to_numpy()) * 100
    return df.hstack(pl.DataFrame(X_qt, schema=gep_qt))


def get_group_mat(df, cat="odor", olfr_key="top_Olfr"):
    olfr = LabelEncoder().fit_transform(df[olfr_key])
    odor = LabelEncoder().fit_transform(df[cat])
    return np.vstack((olfr, odor))


def filter_OR_source_numba(df_OR, SOURCE_THRESH, col="source"):
    """Find ORs present in at least SOURCE_THRESH cells in each of n_sources.

    Parameters
    ----------
    df_OR : [type]
        [description]
    n_sources : [type]
        [OR must be present in at least this_many sources]
    SOURCE_THRESH : int
        [number of ORs for each OR for each of n_sources]
    col : str, optional
        [column of df_OR to use for grouping], by default "source"
    """
    uq_Olfr = np.unique(df_OR.top_Olfr)
    group_mat = get_group_mat(df_OR, cat=col)
    enough_olfr = (npg.aggregate(group_mat, 1, "sum") >= SOURCE_THRESH).all(1)
    good_ORs = uq_Olfr[enough_olfr]
    has_enough_ORs = df_OR.top_Olfr.isin(good_ORs)
    return good_ORs, has_enough_ORs


def subset_source_OR_indices(or_num, THRESH=4, source=None, seed=None, n_boot=100):
    """Subset equal numbers of cells from each source for each OR n_boot times

    Arguments:
        or_num {[type]} -- [np.array of or_labels to randomly sample from]

    Keyword Arguments:
        SOURCE_THRESH {int} -- [number of cells per source per OR] (default: {4})
        source {[type]} -- [vector describing source for each cell] (default: {None})
        seed {[type]} -- [random seed] (default: {None})
        n_boot {int} -- [number different indices to generate] (default: {100})
    """
    np.random.seed(seed)
    tmp = []
    for n in np.unique(or_num):
        if source is None:
            uq_vals = np.where((or_num == n))[0]
            tmp.append(
                uq_vals[np.random.rand(len(uq_vals), n_boot).argsort(0)[:THRESH, :]]
            )
        else:
            for s in np.unique(source):
                uq_vals = np.where((or_num == n) & (source == s))[0]
                tmp.append(
                    uq_vals[np.random.rand(len(uq_vals), n_boot).argsort(0)[:THRESH, :]]
                )
    indices = np.vstack(tmp)
    return indices


def melt_diag(df_corr, k=0, name="value", col_base="env"):
    """Melt upper triangle of symmetric dataframe

    Args:
        df_corr (pd.DataFrame): input dataframe (e.g. )
        k (int, optional): Defaults to 0.

    Returns:
        pd.DataFrame: [description]
    """
    df_corr = df_corr.copy()
    df_corr_nan = df_corr.where(np.triu(np.ones(df_corr.shape), k).astype(bool))
    # make unique names for reset_index()
    if df_corr_nan.index.name == df_corr_nan.columns.name:
        if df_corr_nan.index.name is None:
            df_corr_nan.index.name = col_base + "_1"
            df_corr_nan.columns.name = col_base + "_2"
        else:
            df_corr_nan.index.name = df_corr_nan.index.name + "_1"
            df_corr_nan.columns.name = df_corr_nan.columns.name + "_2"
    df_corr_stack = df_corr_nan.stack()
    df_corr_stack.name = name
    return df_corr_stack.reset_index()
