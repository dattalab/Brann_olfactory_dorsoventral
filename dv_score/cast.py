import numpy as np
import pandas as pd
import polars as pl
from tqdm import trange

from dv_score import io, util
from dv_score.constants import GEP3_qt

WT = "C57BL_6NJ"
KO = "CAST_EiJ"
GENOS = (WT, KO)
C57 = "C57BL/6NJ"
CAST = "CAST/EiJ"
STRAINS = [C57, CAST]
strain_map = dict(zip(GENOS, STRAINS))
DELTA_THRESH = 0.5
F1 = "F$_1$"
F0 = "F$_0$"
F1_THRESH = 2  # at least N cells per OR per strain
GEP_ABS = [g + "_abs" for g in GEP3_qt]
KEYS = ["contig", "chrom", "start", "stop", "ref", "alts"]


def read_update(df, strain_col="strain"):
    df[strain_col] = df[strain_col].map(strain_map)
    return df


def load_files(f1_only=True):
    folders = io.get_data_folders()
    if f1_only:
        df = pd.read_parquet(folders.zenodo / "CAST_C57_F1_df_OR_gep_snp.parquet")
    else:
        df = pd.read_parquet(folders.zenodo / "CAST_C57_F0_F1_df_OR_gep_snp.parquet")

    return read_update(df)


def get_delta(_df, cols=GEP3_qt, l_cols=["top_Olfr", "strain"]):
    sample_delta = (
        _df.group_by(l_cols)
        .mean()
        .sort(l_cols)
        .group_by("top_Olfr")
        .agg([pl.col(col).diff().last() for col in cols])
    )
    return sample_delta.join(
        sample_delta.with_columns([pl.col(col).abs() for col in cols]),
        on="top_Olfr",
        suffix="_abs",
    )


# shuffle strain labels within each OR
def run_shuffle(i, df_pl, cols=GEP3_qt):
    df_subset = (
        df_pl.group_by("top_Olfr")
        .agg([pl.col(col) for col in cols] + [pl.col("strain").shuffle(seed=i)])
        .explode(cols + ["strain"])
    )
    this_delta = get_delta(df_subset, cols=cols)
    return this_delta.with_columns(i=i)


def find_sig_subtypes(df_OR_gep, pct=99):
    folders = io.get_data_folders()
    sig_fn = folders.processed / f"CAST_C57_F1_sig_ORs_pct_{pct}.parquet"
    if sig_fn.exists():
        is_sig = pd.read_parquet(sig_fn)
    else:
        print(f"Calculating sig ORs for each GEP with {pct} % threshold.")
        df_f1 = df_OR_gep[df_OR_gep.abs_strain_delta > DELTA_THRESH].copy()
        _, has_enough = util.filter_OR_source_numba(df_f1, F1_THRESH, col="strain")
        df_pl = pl.DataFrame(df_f1[has_enough])
        all_shuffs = [run_shuffle(i, df_pl) for i in trange(1000, ncols=80)]
        df_shuff = pl.concat(all_shuffs)
        obs_delta = get_delta(df_pl)
        threshes = np.percentile(df_shuff[GEP_ABS], pct, axis=0)
        df_obs = obs_delta.to_pandas().set_index("top_Olfr")
        is_sig = df_obs[GEP_ABS] > threshes
        is_sig.to_parquet(sig_fn)
    return is_sig


def parse_record(rec, anno_col="CSQ", keys=KEYS):
    # only have homozygous variants so can keep first
    gts = [v["GT"] for v in rec.samples.values()]
    vals = [getattr(rec, k) for k in keys]
    return vals + gts + [rec.info[anno_col]]


def make_anno(s):
    return dict([tuple(ss.split(' "')) for ss in s.split('"; ')])
