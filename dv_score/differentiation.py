import polars as pl

from dv_score import io


def load_pseudotime():
    df_obs = io.load_development("cell_obs", cols=["pseudotime_rank"])
    pseudotime = df_obs["pseudotime_rank"]
    x_pseudotime = pseudotime.values
    return pseudotime, x_pseudotime


def med_pseudotime(n_imm_cells=50):
    _, x_pseudotime = load_pseudotime()
    df_olfr_mat = io.load_development("olfr_mat", polars=True).with_columns(
        pseudotime=x_pseudotime
    )
    # only keep cells with pseudotime > 0.4 expressing class II ORs
    # only look at class II ORs with at least `n_imm_cells` cells
    df_olfr_eval = df_olfr_mat.filter(pl.col("pseudotime") > 0.4)
    c2_olfr = set(io.load_olfr_info().query("OR_class == 'class2'").index)
    df_gep_mean = io.load_gep_mean(to_pandas=True)
    c2_keep = c2_olfr.intersection(df_gep_mean.index)
    c2_eval = [c for c in c2_keep if (df_olfr_eval[c] > 0).sum() > n_imm_cells]

    # could also look at percentiles for pseudotime by taking np.percentile(x[x>0], q) for each olfr
    df_pseudo = (df_olfr_eval[c2_eval] > 0).cast(int) * df_olfr_eval["pseudotime"]
    # median pseudotime value of cells expressing OR
    df_plot = (
        df_pseudo.select(pl.col(c2_eval).replace(0, None).median())
        .to_pandas()
        .T.rename(columns={0: "med"})
        .join(df_gep_mean["DV"], how="inner")
    )
    return df_plot
