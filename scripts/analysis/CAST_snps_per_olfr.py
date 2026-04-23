import numpy as np
import pandas as pd
import polars as pl
import pyranges
import pysam
from tqdm import tqdm

from dv_score import io
from dv_score.cast import make_anno, parse_record

bm = io.load_olfr_info()
ens_names = set(bm.ensembl_gene_id)
# see `extract_CAST_olfr_SNPs.sh`
vcf = pysam.VariantFile("mgp_REL2021_snps.olfr.CAST.vcf.gz")
keys = ["contig", "chrom", "start", "stop", "ref", "alts"]
SAMPLES = list(vcf.header.samples)
col_names = keys + SAMPLES + ["anno"]
print(SAMPLES)
rec = next(iter(vcf))
C57 = SAMPLES[0]
CAST = SAMPLES[1]
WT_SNP = (0, 0)
MUT_SNP = (1, 1)


vcf = pysam.VariantFile(
    "snpEff_annotated/mgp_REL2021_snps.olfr.CAST.clean.ann.vcf.gz"
)
all_recs = []
for rec in tqdm(vcf):
    all_recs.append(parse_record(rec, anno_col="ANN"))

df_vcf_anno = pd.DataFrame(all_recs, columns=col_names).convert_dtypes()
is_anno_homo = df_vcf_anno[CAST].isin([(0,0), (1,1), (2,2)])
df_vcf_anno = df_vcf_anno[is_anno_homo].copy()
df_vcf_anno["has_missense"] = df_vcf_anno["anno"].apply(
    lambda l: any(["missense" in i for i in l])
)

df_olfr_gtf = (
    pl.scan_csv(
        "Mus_musculus.GRCm39.105.with_manual_taar_utr_sorted.gtf",
        separator="\t",
        comment_prefix="#",
        ignore_errors=True,
        has_header=False,
        infer_schema=False,
        new_columns=[
            "seqname",
            "source",
            "type",
            "start",
            "end",
            "v1",
            "strand",
            "v2",
            "desc",
        ],
    )
    # keep only CDS annotations that are in OR genes
    .filter(pl.col("type").is_in(("CDS", "exon", "three_prime_utr", "five_prime_utr")))
    .filter(
        pl.col("desc")
        .str.split(";")
        .list.get(0)
        .str.split('"')
        .list.get(1)
        .is_in(ens_names)
    )
    .with_columns(pl.col(("start", "end")).cast(int))
    .with_columns(size=pl.col("end") - pl.col("start"))
    .unique(("seqname", "type", "start", "end"))
    .sort("seqname", "start")
).collect()


bed_cols = ["Chromosome", "Start", "End"]
df_vcf_bed = df_vcf_anno[["chrom", "start", "stop", "has_missense"]]
df_vcf_bed.columns = bed_cols + ["has_missense"]
df_vc_pr = pyranges.PyRanges(df_vcf_bed)

df_cds_gtf = df_olfr_gtf.to_pandas()
df_cds_anno = df_cds_gtf.join(
    pd.DataFrame(df_cds_gtf["desc"].apply(make_anno).to_dict()).T
)
df_cds_merge = df_cds_anno.merge(
    bm.drop(columns="external_gene_name").reset_index(),
    left_on=["gene_id"],
    right_on="ensembl_gene_id",
)
print((df_cds_merge["external_gene_name"] == df_cds_merge["gene_name"]).mean())

new_bed_cols = bed_cols + ["top_Olfr"]
df_cds_bed = df_cds_merge[["seqname", "start", "end", "external_gene_name"]]
df_cds_bed.columns = new_bed_cols
df_cds_pr = pyranges.PyRanges(df_cds_bed)
df_cds_uq = df_cds_pr

# number of CDS snps per OR
is_cds = df_pr_inter = df_cds_gtf.type == "CDS"
df_pr_inter = df_cds_uq[is_cds].intersect(df_vc_pr)
n_cds_snps = df_pr_inter.df.groupby("top_Olfr")["Start"].nunique()

# number of missense variants per OR
miss_inter = df_cds_uq[is_cds].intersect(df_vc_pr[df_vc_pr.has_missense])
n_miss_per_olfr = miss_inter.df.groupby("top_Olfr")["Start"].nunique()

overlaps = df_cds_uq.count_overlaps(df_vc_pr)
assert len(overlaps) == df_cds_merge.shape[0]

df_cds_overlap = df_cds_merge.merge(
    overlaps.df.drop_duplicates(),
    left_on=["seqname", "start", "end", "external_gene_name"],
    right_on=new_bed_cols,
    how="left",
)

df_n_per = df_cds_overlap.groupby(["type", "top_Olfr"], as_index=False)[
    ["size", "NumberOverlaps"]
].sum()
df_n_per["snp_per_kb"] = df_n_per["NumberOverlaps"] / df_n_per["size"] * 1e3

# SNPs across gene level
df_olfr_gene = (
    pl.scan_csv(
        "Mus_musculus.GRCm39.105.with_manual_taar_utr_sorted.gtf",
        separator="\t",
        comment_prefix="#",
        ignore_errors=True,
        has_header=False,
        infer_schema=False,
        new_columns=[
            "seqname",
            "source",
            "type",
            "start",
            "end",
            "v1",
            "strand",
            "v2",
            "desc",
        ],
    )
    # keep only CDS annotations that are in OR genes
    .filter(pl.col("type") == "gene")
    .with_columns(
        ensembl_gene_id=pl.col("desc")
        .str.split(";")
        .list.get(0)
        .str.split('"')
        .list.get(1)
    )
    .filter(pl.col("ensembl_gene_id").is_in(ens_names))
).collect()
df_bm_all = bm.merge(df_olfr_gene.to_pandas(), on="ensembl_gene_id", how="left")

# get strand
# if opposite strand, flip coordinates
se = bm[["start_position", "end_position"]]
inds = np.array([[1, 0], [0, 1]])
to_flip = ((bm.Strand + 1) / 2).astype(int).values
df_flip = pd.DataFrame(
    se.values[np.arange(se.shape[0])[:, None], inds[to_flip]],
    index=bm.index,
    columns=["strand_start", "strand_end"],
)
bm_with_strand = bm.join(df_flip)

# number of
is_pos_strand = bm_with_strand.Strand == 1
bm_pos = bm_with_strand[is_pos_strand].copy()
bm_neg = bm_with_strand[~is_pos_strand].copy()
offset = 1000
bm_pos["up_start"] = bm_pos["start_position"] - offset
# gene_end is 5' side for negative strand
bm_neg["down_end"] = bm_neg["end_position"] + offset
bm1 = bm_pos[["chromosome_name", "up_start", "start_position"]].copy()
bm1.columns = bed_cols
bm2 = bm_neg[["chromosome_name", "end_position", "down_end"]].copy()
bm2.columns = bed_cols

df_utr_bed = pyranges.PyRanges(
    pd.concat([bm1, bm2]).sort_values(["Chromosome", "Start"]).reset_index()
)

df_upstream_overlap = df_utr_bed.count_overlaps(df_vc_pr)
to_join = (
    df_upstream_overlap.df[["external_gene_name", "NumberOverlaps"]]
    .copy()
    .rename(columns={"external_gene_name": "top_Olfr"})
)
to_join["type"] = "Promoter"
to_join["size"] = offset
to_join["snp_per_kb"] = to_join["NumberOverlaps"] / to_join["size"] * 1e3

df_n_per_all = pd.concat([df_n_per, to_join[df_n_per.columns]])
