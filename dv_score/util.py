from collections import namedtuple
from pathlib import Path

import pandas as pd
import polars as pl


def get_data_path():
    # data folder in base directory
    return Path(__file__).parents[1] / "data"


def get_data_folders(data_path=None):
    if data_path is None:
        data_path = get_data_path()
    fields = ["processed", "tables"]
    dir = namedtuple("dir", fields)
    data_fold = dir(*[data_path / f for f in fields])
    return data_fold


def load_olfr_info():
    data_fold = get_data_folders()
    return pd.read_csv(
        data_fold.tables / "Olfr_biomart_mm39_Ensembl_105.csv", index_col=0
    )


def load_df_dv(filter=True, vc_thresh=150):
    folders = get_data_folders()
    df_mature = pl.read_parquet(folders.processed / "DV_score_all_mature_OSNs.parquet")
    if filter:
        vc = df_mature["top_Olfr"].value_counts()
        uq_olfr = (
            vc.filter(pl.col("count") >= vc_thresh).select(pl.col("top_Olfr")).unique()
        )
        df_mature = df_mature.filter(pl.col("top_Olfr").is_in(uq_olfr))
    return df_mature
