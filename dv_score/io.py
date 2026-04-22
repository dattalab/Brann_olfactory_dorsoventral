import tempfile
import urllib.request
from collections import namedtuple
from pathlib import Path

import pandas as pd
import polars as pl

from dv_score import constants, util


def get_data_path():
    """Return the path to the data directory in the base project directory."""
    return Path(__file__).parents[1] / "data"


def get_data_folders(data_path=None):
    """Return a named tuple of key data subdirectory paths.

    Parameters
    ----------
    data_path : Path or str, optional
        Root data directory. Defaults to the project data directory.

    Returns
    -------
    dir : namedtuple
        Named tuple with fields ``processed``, ``zenodo``, and ``tables``.
    """
    if data_path is None:
        data_path = get_data_path()
    dir = namedtuple("dir", ["processed", "zenodo", "tables"])
    folds = ["processed", "processed/from_zenodo", "tables"]
    data_fold = dir(*[data_path / f for f in folds])
    return data_fold


def load_olfr_info():
    """Load olfactory receptor gene info from BioMart (mm39, Ensembl 105).

    Returns
    -------
    pd.DataFrame
        DataFrame with Olfr gene metadata, exported from Ensembl Biomart
    """
    data_fold = get_data_folders()
    return pd.read_csv(
        data_fold.tables / "Olfr_biomart_mm39_Ensembl_105.csv", index_col=0
    )


def load_df_dv(apply_filter=True, vc_thresh=150, add_qt=True, qt_geps=["DV"]):
    """Load DV score data for mature OSNs.

    Parameters
    ----------
    apply_filter : bool, optional
        If True, retain only OR subtypes with at least ``vc_thresh`` cells.
    vc_thresh : int, optional
        Minimum cell count per OR subtype when filtering. Default 150.
    add_qt : bool, optional
        If True, append quantile-transformed GEP columns via ``util.pl_qt``.
    qt_geps : list of str, optional
        GEP column names to quantile-transform. Default ``["DV"]``.

    Returns
    -------
    pl.DataFrame
        Polars DataFrame of mature OSN DV scores, optionally filtered and
        quantile-transformed.
    """
    folders = get_data_folders()
    df_mature = pl.read_parquet(folders.processed / "DV_score_all_mature_OSNs.parquet")
    if apply_filter:
        vc = df_mature["top_Olfr"].value_counts()
        uq_olfr = (
            vc.filter(pl.col("count") >= vc_thresh).select(pl.col("top_Olfr")).unique()
        )
        df_mature = df_mature.join(uq_olfr, on="top_Olfr", how="inner")
        # df_mature = df_mature.filter(pl.col("top_Olfr").is_in(uq_olfr))
    if add_qt:
        df_mature = util.pl_qt(df_mature, geps=qt_geps)
    return df_mature


def load_gep_mean(folders=None, vc_thresh=150, to_pandas=False):
    """Load mean GEP usage per OR subtype from Zenodo data.

    Parameters
    ----------
    folders : namedtuple, optional
        Data folder paths. Defaults to ``get_data_folders()``.
    vc_thresh : int, optional
        Minimum cell count per OR subtype to include. Default 150.
    to_pandas : bool, optional
        If True, return a pandas DataFrame indexed by ``top_Olfr``.

    Returns
    -------
    pl.DataFrame or pd.DataFrame
        Mean GEP usage per OR subtype, filtered by number of OSNs per subtype.
    """
    if folders is None:
        folders = get_data_folders()
    df_gep_mean = (
        pl.read_parquet(folders.zenodo / "OSN_subtype_mean_GEP_usage.parquet")
        .filter(pl.col("n_cells") >= vc_thresh)
        .drop(["n_cells"])
    )
    if to_pandas:
        df_gep_mean = df_gep_mean.to_pandas().set_index("top_Olfr")
    return df_gep_mean


def load_hvg_mean():
    """Load mean expression of highly variable genes per OR subtype.

    Returns
    -------
    pd.DataFrame
        DataFrame of mean HVG expression per OR subtype.
    """
    folders = get_data_folders()
    return pd.read_csv(
        folders.zenodo / "OSNDev.Data_1_hvg_mean_per_OR.csv.gz", index_col=0
    )


def load_development(key, cols=None, polars=False):
    """Load INP5 OSN development data by key.

    Parameters
    ----------
    key : str
        Dataset to load. One of ``"cell_obs"``, ``"olfr"``, or ``"olfr_mat"``.
    cols : list of str, optional
        Subset of columns to load. If None, all columns are returned.
    polars : bool, optional
        If True, return a Polars DataFrame; otherwise return pandas.

    Returns
    -------
    pl.DataFrame or pd.DataFrame
        Development data for the requested key.
    """
    f_map = {
        "cell_obs": "INP5_df_obs.parquet",
        "olfr": "INP5_df_olfr.parquet",
        "olfr_mat": "INP5_df_olfr_mat.parquet",
    }
    if key not in f_map:
        raise ValueError(f"Key {key} must be one of {f_map.keys()}")
    folders = get_data_folders()
    if polars:
        return pl.read_parquet(folders.zenodo / f_map[key], columns=cols)
    else:
        return pd.read_parquet(folders.zenodo / f_map[key], columns=cols)


def load_clones():
    """Load OSN lineage tracing (clonal) data.

    Returns
    -------
    pd.DataFrame
        DataFrame of OSN lineage tracing data.
    """
    folders = get_data_folders()
    return pd.read_parquet(folders.zenodo / "OSN_lineage_tracing.parquet")


def load_ra(rar=False):
    """Load retinoic acid (RA) or RAR inhibition DV score data."""
    folders = get_data_folders()
    if rar:
        return pd.read_parquet(folders.zenodo / "RAR_inh_df_wt_dv.parquet")
    else:
        return pd.read_parquet(folders.zenodo / "RA_inh_df_wt_dv.parquet")


def load_glomeruli(return_pandas=True):
    """Load olfactory bulb Visium glomeruli positions joined with GEP mean data.

    Parameters
    ----------
    return_pandas : bool, optional
        If True, return a pandas DataFrame; otherwise return a Polars DataFrame.

    Returns
    -------
    pd.DataFrame or pl.DataFrame
        Glomeruli spatial positions joined with mean GEP usage per OR subtype.
    """
    folders = get_data_folders()
    df_gep_qt_mean = (
        load_df_dv(qt_geps=constants.GEP3)
        .group_by("top_Olfr")
        .agg(pl.col(constants.GEP3_qt + constants.GEP3).mean())
    )
    # df_gep_mean = load_gep_mean(folders)
    df_geom_dv = pl.read_parquet(
        folders.zenodo / "Olfactory_bulb_visium_glomeruli_positions.parquet"
    ).join(df_gep_qt_mean, on="top_Olfr", how="inner")
    if return_pandas:
        df_geom_dv = df_geom_dv.to_pandas()
    return df_geom_dv


def load_other_spatial(key):
    """Load DV score comparisons against other spatial datasets.

    Parameters
    ----------
    key : str
        Dataset to load. One of ``"tan"`` (zone index), ``"saraiva"``
        (3D index), or ``"bintu"`` (MERFISH).

    Returns
    -------
    pd.DataFrame
        DataFrame for the requested spatial comparison dataset.
    """
    f_map = {
        "tan": "DV_score_3M_vs_tan_inferred_zone_index.parquet",
        "saraiva": "DV_score_3M_vs_saraiva_3D_index.parquet",
        "bintu": "MERFISH_KNN_regress_df_pred_all.parquet",
    }
    if key not in f_map:
        raise ValueError(f"Key {key} must be one of {f_map.keys()}")
    folders = get_data_folders()
    return pd.read_parquet(folders.zenodo / f_map[key])


def load_stereo_seq():
    """Load Stereo-seq DV KNN classification and binned spatial data.

    Returns
    -------
    df : pd.DataFrame
        Cell-level DV KNN classification results for two Stereo-seq sections.
    df_sum_rot : pd.DataFrame
        Binned mean expression data for one section, aligned and with isolated
        spots removed.
    """
    folders = get_data_folders()
    df = pd.read_parquet(
        folders.zenodo / "C06128C3_C06112G5_DV_KNN_classification.parquet",
        engine="fastparquet",
    )
    # load binned mean
    df_sum_rot = pd.read_parquet(folders.zenodo / "C06128C3_df_sum_rot_aligned.parquet")
    df_sum_rot = df_sum_rot[~df_sum_rot.isolated]
    return df, df_sum_rot


def load_merfish():
    """Load MERFISH coronal section DV score predictions and spatial alignment data.

    Cells are filtered to those present in both the prediction and alignment
    DataFrames.

    Returns
    -------
    df_pred_keep : pd.DataFrame
        Predicted DV scores for cells that appear in the alignment data.
    df_align : pd.DataFrame
        Cell annotations with polygon-based spatial coordinates and DV alignment
        scores, including an added ``pwl_rescale_invert`` column.
    """
    folders = get_data_folders()
    df_pred = pd.read_parquet(
        folders.zenodo / "Coronal2_DBSCAN_cells_predicted_DV_score_global_knn.parquet"
    )
    df_align = pd.read_parquet(
        folders.zenodo
        / "Coronal2_cells_annotated_per_polygon_rank_distance_common_rank_pwl_regression_drop_dups.parquet"
    )
    df_align["pwl_rescale_invert"] = 1 - df_align["pwl_rescale"]

    cell_keep = df_pred.index.isin(df_align.olfr_cell_labels)
    df_pred_keep = df_pred[cell_keep]

    return df_pred_keep, df_align


def load_merfish_AP():
    """Load MERFISH apical-basal distance and AP score data."""
    folders = get_data_folders()
    _df = pd.read_parquet(
        folders.zenodo / "MERFISH_apical_basal_distance_AP_score.parquet"
    )
    return _df[_df.close_dist]


def load_bintu_index():
    url = (
        "https://raw.githubusercontent.com/"
        "BogdanBintu/MERFISH_Analysis_Olfactory_Receptors/master/"
        "Processed_data/Table%20S3.xlsx"
    )

    with tempfile.NamedTemporaryFile(suffix=".xlsx") as tmp:
        urllib.request.urlretrieve(url, tmp.name)
        df_index = pd.read_excel(tmp.name, comment="#")
        df_index.columns = df_index.columns.map(
            lambda l: l.replace(" ", "_").lower().replace("-", "_")
        )
        df_index = df_index.set_index("receptor_name")
        df_index = df_index[~df_index.central_to_peripheral_index.isnull()]
    return df_index
