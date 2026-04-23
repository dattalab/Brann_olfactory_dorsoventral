import numpy as np
import pandas as pd
import pyranges
import pysam
from joblib import Parallel, delayed
from tqdm import tqdm

from dv_score.cast import C57, CAST, KEYS, load_files, parse_record
from dv_score.util import get_cores

N_READ_THRESH = 1
READ_FRAC_THRESH = 0.8
DELTA_THRESH = 0.5


def get_snps_in_olfr(in_tup, tags=["RG", "CB", "UB"], attrs=["reference_name"]):
    olfr, row, df_this_olfr = in_tup
    # account for Olfr290 with multiple entries in bed file
    if len(row.shape) > 1:
        chrom = row.Chromosome[0]
        start = row.Start.min()
        stop = row.End.max()
    else:
        chrom, start, stop = row
    # assert (df_this_olfr.chrom == chrom).all()
    # assert (df_this_olfr.Name == olfr).all()
    vcf_pos = set(df_this_olfr.start)

    olfr_bam = "olfr_bam_subset/F1_mm39_olfr.bam"
    bam = pysam.AlignmentFile(olfr_bam)

    outs = []
    pile = bam.fetch(chrom, start, stop, multiple_iterators=True)
    for read in pile:
        pos_overlap = set(read.get_reference_positions()) & vcf_pos
        if len(pos_overlap) > 0:
            r_tags = [read.get_tag(t) for t in tags] + [getattr(read, a) for a in attrs]
            seq = read.query_sequence
            for s, pos in read.get_aligned_pairs(matches_only=True):
                if pos in pos_overlap:
                    outs.append([seq[s], read.query_qualities[s], pos] + r_tags)

    df_out = pd.DataFrame(outs, columns=["seq", "qual", "start"] + tags + attrs)
    df_out["top_Olfr"] = olfr
    return df_out


df_pr = pyranges.read_bed("mm39_olfr.bed")
df_pr_df = df_pr.df.set_index("Name")
all_recs = []
# get all CAST snps in OR regions
# see `extract_CAST_olfr_SNPs.sh`
vcf = pysam.VariantFile("mgp_REL2021_snps.olfr.CAST.vcf.gz")
SAMPLES = list(vcf.header.samples)
for rec in tqdm(vcf):
    all_recs.append(parse_record(rec))
col_names = KEYS + SAMPLES + ["anno"]
df_vcf = pd.DataFrame(all_recs, columns=col_names).convert_dtypes()
is_homo = (df_vcf[CAST] == (2, 2)) | (df_vcf[CAST] == (1, 1))
df_vcf_homo = df_vcf[is_homo].copy()
df_vcf_homo["alt"] = df_vcf_homo.apply(lambda l: l.alts[l[CAST][0] - 1], axis=1)

# vcf homozygous snps
bed_cols = ["Chromosome", "Start", "End"]
df_vcf_bed = df_vcf_homo[["chrom", "start", "stop"]]
df_vcf_bed.columns = bed_cols
df_vc_pr = pyranges.PyRanges(df_vcf_bed)

# intersect OR bed file with homo SNP vcf
df_pr_name = df_pr.intersect(df_vc_pr)
# get gene name for each SNP
df_vcf_merge = df_vcf_homo.merge(
    df_pr_name.df, left_on=["chrom", "start", "stop"], right_on=bed_cols, how="left"
)
# some SNPs correspond to multiple genes
# try to avoid duplicate names for ORs in the same location
df_OR_gep = load_files()
olfr_found = set(df_OR_gep.top_Olfr)
df_vcf_merge_found = df_vcf_merge[df_vcf_merge.Name.isin(olfr_found)]
nu = df_vcf_merge_found.groupby(["chrom", "start"]).Name.nunique()
# keep SNPs only matching one OR
good_pos = nu[nu == 1].index
df_vcf_keep = (
    df_vcf_merge_found.set_index(["chrom", "start"]).loc[good_pos].reset_index()
)
assert pd.Index(df_vcf_keep.Name.unique()).isin(df_pr_df.index).all()
in_tups = []
for olfr, df_this_vcf in df_vcf_keep.groupby("Name"):
    in_tups.append((olfr, df_pr_df.loc[olfr], df_this_vcf))

# run for each OR
cores = get_cores()
results = Parallel(n_jobs=cores, verbose=True)(
    delayed(get_snps_in_olfr)(i) for i in in_tups
)
df_seqs = pd.concat(results)
print(df_seqs.shape)
df_seq_merge = df_seqs.merge(
    df_vcf_keep[["chrom", "start", "ref", "alt"]],
    left_on=["reference_name", "start"],
    right_on=["chrom", "start"],
).drop(columns=["reference_name"])
df_seq_merge["cell_barcode"] = (
    df_seq_merge["RG"] + "_" + df_seq_merge["CB"].apply(lambda l: l.split("-")[0])
)
df_seq_merge["is_ref"] = df_seq_merge["seq"] == df_seq_merge["ref"]
df_seq_merge["is_alt"] = df_seq_merge["seq"] == df_seq_merge["alt"]
df_seq_merge.to_parquet("F1_mm39_olfr_snp_seq.parquet")

# filter for only reads in chosen OR
# keep only reads in chosen OR
df_seq_inner = df_seq_merge.merge(
    df_OR_gep.reset_index()[["index", "top_Olfr", "OR_counts"]],
    left_on=["cell_barcode", "top_Olfr"],
    right_on=["index", "top_Olfr"],
    how="inner"
)
# only 1 OR per cell barcode
assert (df_seq_inner.groupby("cell_barcode").top_Olfr.nunique() == 1).all()

is_cols = ["is_ref", "is_alt"]
qual_thresh = 30
df_qual = df_seq_inner[df_seq_inner.qual >= qual_thresh]
# n_snps_per_olfr = df_qual.groupby("top_Olfr", as_index=False).start.nunique()

n_bc_reads = df_qual.cell_barcode.value_counts().rename("n_reads")
# number of fraction of SNPs for each strain for each OR
ds = []
for agg_func in [np.sum, np.mean]:
    d = df_qual.groupby(["cell_barcode"])[is_cols].agg(agg_func)
    d["either"] = d[is_cols].sum(1)
    d["max"] = d[is_cols].max(1)
    ds.append(d)
d_summary = ds[0].join(ds[1], rsuffix="_frac").join(n_bc_reads)
d_summary["strain_delta"] = d_summary["is_ref_frac"] - d_summary["is_alt_frac"]
d_summary["abs_strain_delta"] = d_summary["strain_delta"].abs()
d_summary["ref_higher"] = d_summary.strain_delta > 0
d_summary["strain_max"] = d_summary["ref_higher"].map({False: CAST, True: C57})

d_keep = d_summary[
    (d_summary.n_reads >= N_READ_THRESH) & (d_summary.max_frac >= READ_FRAC_THRESH)
]
df_OR_gep = df_OR_gep.join(d_keep, how="left", suffix="_snp")
