
import itertools
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

from dv_score import util


def flatten(input_list):
    return [item for sublist in input_list for item in sublist]


def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def hamming(bc1, bc2):
    return np.sum([x1 != x2 for x1, x2 in zip(bc1, bc2)])


def parallel_process(array, function, n_jobs=16, front_num=0, process=True, **kwargs):
    """
        A parallel version of the map function with a progress bar. 

        Args:
            array (array-like): An array to iterate over.
            function (function): A python function to apply to the elements of array
            n_jobs (int, default=16): The number of cores to use
            use_kwargs (boolean, default=False): Whether to consider the elements of array as dictionaries of 
                keyword arguments to function 
            front_num (int, default=3): The number of iterations to run serially before kicking off the parallel job. 
                Useful for catching bugs
        see http://danshiebler.com/2016-09-14-parallel-progress-bar/
        Returns:
            [function(array[0]), function(array[1]), ...]
    """
    if process:
        processor = ProcessPoolExecutor
    else:
        processor = ThreadPoolExecutor

    # We run the first few iterations serially to catch bugs
    # tq_use, tr_use = get_tqdm()
    if front_num > 0:
        front = [function(a, **kwargs) for a in array[:front_num]]
    else:
        front = []
    # If we set n_jobs to 1, just run a list comprehension. This is useful for benchmarking and debugging.
    if n_jobs == 1:
        return front + [function(a, **kwargs) for a in tqdm(array[front_num:], leave=True)]
    # Assemble the workers
    with processor(max_workers=n_jobs) as pool:
        # Pass the elements of array into function
        futures = [pool.submit(function, a, **kwargs) for a in array[front_num:]]
        tqdm_kwargs = {"total": len(futures), "unit": "it", "unit_scale": True, "leave": True}
        # Print out the progress as tasks complete
        for _ in tqdm(as_completed(futures), **tqdm_kwargs):
            pass
    out = []
    # Get the results from the futures.
    for _, future in enumerate(futures):
        try:
            out.append(future.result())
        except Exception as e:
            out.append(e)
    return front + out


def get_cell_barcodes(tup_in, umi_thresh=4, bc_thresh=3, bc_frac=0.75):
    cell, _df = tup_in
    umi_c = _df.UB_collapse.value_counts()
    umi_eval = umi_c[umi_c >= umi_thresh]
    bc_found = {}
    for umi in umi_eval.index:
        this_umi = _df[_df.UB_collapse == umi]
        umi_bc_c = this_umi.bc2_collapse.value_counts()
        top_counts = umi_bc_c.max()
        try:
            if top_counts >= bc_thresh and top_counts >= (bc_frac * umi_bc_c.sum()):
                top_bc2 = umi_bc_c.idxmax()
                top_bc1 = (
                    this_umi[this_umi.bc2_collapse == top_bc2]
                    .bc1_correct.value_counts()
                    .index[0]
                )
                bc_found[(cell, umi)] = (top_bc1, top_bc2, top_counts, umi_c.loc[umi])
        # account for b1_correct sometimes being empty when BC1 isn't found
        except IndexError as err:
            pass
    return bc_found


def cell_barcodes_frac(tup_in, umi_thresh=4, bc_thresh=3, bc_frac=0.75):
    cell, _df = tup_in
    umi_c = _df.UB_collapse.value_counts()
    umi_eval = umi_c[umi_c >= umi_thresh]
    bc_found = []
    for rank, umi in enumerate(umi_eval.index):
        this_umi = _df[_df.UB_collapse == umi]
        umi_bc_c = this_umi.bc2_collapse.value_counts()
        top_bc2 = umi_bc_c.idxmax()
        top_counts = umi_bc_c.max()
        try:
            top_bc1 = (
                this_umi[this_umi.bc2_collapse == top_bc2].bc1_correct.mode().loc[0]
            )
        except KeyError as err:
            top_bc1 = (
                this_umi[this_umi.bc2_collapse == top_bc2]
                .bc1_correct.astype(str)
                .mode()
                .loc[0]
            )

        sorted_bc_counts = umi_bc_c.values
        above = np.cumsum(sorted_bc_counts[1:]) > top_counts
        if above.any():
            n_above = np.argmax(above)
        else:
            n_above = np.nan
        #     assert sorted_bc_counts.sum() == umi_c.loc[umi]
        if len(sorted_bc_counts) == 1:
            tops = (None, None, None, None, None)
        else:
            tops = (
                n_above,
                top_counts / sorted_bc_counts[1],
                top_counts / sorted_bc_counts.sum(),
                top_counts / sorted_bc_counts[:5].sum(),
                top_counts / sorted_bc_counts[:10].sum(),
            )

        bc_found.append(
            (cell, umi, rank, top_bc1, top_bc2, top_counts, umi_c.loc[umi], *tops)
        )
    return bc_found


MOUSE = "cellecta_round7_GFPpos"
save_data_folder = ""
adata = "" # combined adata object from lentiviral experiments
cell1 = adata.obs_names[adata.obs.orig_ident == MOUSE].map(
    lambda l: l.split("_")[-1]
)
cell_set = set(cell1)

next_folder = "" # nextflow pipeline
cell_folder = next_folder / f"{MOUSE}/cell_ranger/{MOUSE}/outs"
cell_folder2 = next_folder / f"{MOUSE}/"
PRIMER = "CCGACCACCGAACGCAACGCACGCA"
XP = "TGGTGCAGCTGGAGCAG"

primer_len = len(PRIMER)
xp_len = len(XP)
AMP_LEN = 90
bc1_len = 14
bc2_len = 30
constant = "TGGT"
constant_len = len(constant)

# see `extract_barcode and nextflow+_with_clonal_barcode.nf`
samfile = pysam.AlignmentFile(cell_folder2 / "possorted_barcode_with_header.bam", "rb")
reads = []
for read in tqdm(samfile):
    try:
        reads.append((read.seq, read.get_tag("CB"), read.get_tag("UB")))
    except KeyError as err:
        pass

df_reads = pd.DataFrame(reads, columns=["seq", "CB", "UB"])
print(df_reads.shape)
df_reads["cell"] = df_reads.CB.map(lambda l: l.split("-")[0])
print(df_reads.CB.nunique())
df_reads["primer_start"] = df_reads.seq.apply(lambda l: l.startswith(PRIMER))
print(df_reads["primer_start"].mean() * 100)
df_reads["right_end"] = df_reads.seq.apply(lambda l: l[:AMP_LEN].endswith(XP))
df_reads["has_end"] = df_reads.seq.apply(lambda l: XP in l)
print(df_reads.right_end.mean() * 100)
in_cell = df_reads["cell"].isin(cell_set)
df_reads_cell = df_reads[(in_cell) & (df_reads.primer_start)].copy()
df_reads_cell_good = df_reads_cell[df_reads_cell.right_end]
print(df_reads_cell_good.shape)

lens = [0, primer_len, bc1_len, constant_len, bc2_len, xp_len]
split_pos = list(itertools.accumulate(lens))
pairs_use = list(pairwise(split_pos))[1:-1]
print(pairs_use)

# get BC14 and BC30
splits = []
for s in tqdm(df_reads_cell_good.seq):
    splits.append([s[i:j] for i, j in pairs_use])

df_bcs = pd.DataFrame(
    splits, index=df_reads_cell_good.index, columns=["BC1", "mid", "BC2"]
)
df_reads_cell_good_bc = df_reads_cell_good.join(df_bcs)
df_reads_cell_good_bc["good_mid"] = df_reads_cell_good_bc.mid == constant
print(df_reads_cell_good_bc["good_mid"].mean() * 100)

assert df_reads_cell_good_bc.index.equals(df_reads_cell_good.index)
reads_per_cell = df_reads_cell_good.cell.value_counts()
thresh = [1, 2, 3, 5, 10, 50, 100, 1000, 2000, 5_000, 10_000, 100_000]
ns = {t: (reads_per_cell >= t).sum() for t in thresh}
print(ns)


# load cellecta library barcodes
barcode_folder = save_data_folder / "barcodes/Cellecta-SEQ-CloneTracker-XP-10M-BC14-BC30-Barcodes"
df_bc14 = pd.read_excel(
    barcode_folder / "Cellecta-SEQ-CloneTracker-XP-Barcode-Libraries-100_x_BC14.xlsx",
    skiprows=9,
)
df_bc14.columns = ["Barcode", "ID", "sense", "anti", "1M", "10M"]
df_bc14 = df_bc14.set_index('ID')
df_bc30 = pd.read_excel(
    barcode_folder / "Cellecta-SEQ-CloneTracker-XP-Barcode-Libraries-100K_x_BC30.xlsx",
    skiprows=9,
)
df_bc30.columns = ["Barcode", "ID", "sense", "anti"]
df_bc30 = df_bc30.set_index("ID")
bc30_bcs = set(df_bc30.sense)
bc14_bcs = set(df_bc14.sense)


# look at umis with at least 10 reads per UMI
reads_per_umi = df_reads_cell.UB.value_counts()
umi_subset = reads_per_umi.index[(reads_per_umi >= 10)]
print(len(umi_subset))
df_subset = df_reads_cell_good_bc[df_reads_cell_good_bc.UB.isin(umi_subset)]

# correct bc14 and bc30 based on cellecta
uq_bc1 = df_reads_cell_good_bc.BC1.unique()
df_bc1 = pd.DataFrame(
    [u in bc14_bcs for u in uq_bc1], index=uq_bc1, columns=["is_bc14"]
)
df_bc1 = df_bc1.join(df_reads_cell_good_bc.BC1.value_counts().rename("BC1_reads"))

uq_bc2 = df_reads_cell_good_bc.BC2.unique()
# bcs with perfect match to white list
df_bc2 = pd.DataFrame(
    [u in bc30_bcs for u in uq_bc2], index=uq_bc2, columns=["is_bc30"]
)
df_bc2 = df_bc2.join(df_reads_cell_good_bc.BC2.value_counts().rename("BC2_reads"))
found_bcs = df_bc1.index[df_bc1.is_bc14]
other_bcs = df_bc1.index[~df_bc1.is_bc14]
# hamming distance between two cellecta barcodes
hams = [hamming(i, j) for (i, j) in itertools.combinations(found_bcs, 2)]

bc14_thresh = 3
bc_mapping = {}
n_two = 0
for bc in tqdm(other_bcs):
    n_map = 0
    new_bc = None
    for test_bc in bc14_bcs:
        if hamming(bc, test_bc) <= bc14_thresh:
            n_map += 1
            new_bc = test_bc
        if n_map == 1:
            bc_mapping[bc] = new_bc
        if n_map > 1:
            n_two += 1
for bc in found_bcs:
    bc_mapping[bc] = bc
df_bc1["bc1_correct"] = df_bc1.index.map(bc_mapping)
df_bc1["correct_found"] = df_bc1["bc1_correct"].isin(bc14_bcs)
print(n_two)
print(df_bc1.correct_found.mean())

found_bc_30 = set(df_bc2.index[df_bc2.is_bc30])
other_bc_30 = df_bc2.index[~df_bc2.is_bc30]

hams2 = [hamming(i, j) for (i, j) in itertools.combinations(found_bc_30, 2)]

bc30_thresh = 5


def check_hamming(bc, found=found_bc_30, hamming_thresh=bc30_thresh):
    n_map = 0
    for test_bc in found_bc_30:
        if hamming(bc, test_bc) <= bc30_thresh:
            n_map += 1
            new_bc = test_bc
    if n_map == 1:
        return (bc, new_bc)


cores = util.get_cores()
bc_maps = parallel_process(other_bc_30, check_hamming, n_jobs=cores, front_num=0)

bc_30_mapping = dict([b for b in bc_maps if b is not None])
for bc in found_bc_30:
    bc_30_mapping[bc] = bc
df_bc2["bc2_correct"] = df_bc2.index.map(bc_30_mapping)
df_bc2["correct_found"] = df_bc2["bc2_correct"].isin(bc30_bcs)
df_bc2["correct_found"].mean() * 100

# correct BCs where they have lots of reads but don't perfectly match any cellecta bc
# likely pcr errors early on?
no_match = other_bc_30.difference(bc_30_mapping.keys())
print(len(no_match))
# sort ranked no_match
new_match = {}
df_no_match = df_bc2.loc[no_match].sort_values(['BC2_reads'], ascending=False)
for i, bc in tqdm(enumerate(df_no_match.index)):
    # for each barcode go from n+1:end for 
    if bc in new_match:
        continue
    test_bcs = df_no_match.index[i+1:].difference(new_match.keys())
    for test in test_bcs:
        if hamming(bc, test) <= bc30_thresh:
            new_match[test] = bc
# add back original ones
for bc in df_no_match.index:
    if bc not in new_match:
        new_match[bc] = bc
df_no_match["new_bc"] = df_no_match.index.map(new_match)

bc30_either_map = {**bc_30_mapping, **new_match}
df_bc2["bc2_collapse"] = df_bc2.index.map(bc30_either_map)

assert (
    df_bc2[df_bc2.correct_found].bc2_correct
    == df_bc2[df_bc2.correct_found].bc2_collapse
).all()

# keep reads with BC1 and BC2
df_bc1_keep = df_bc1[df_bc1.correct_found]
bc1_mapping = df_bc1_keep.bc1_correct.to_dict()
df_reads_cell_good_bc["bc1_correct"] = df_reads_cell_good_bc.BC1.map(bc1_mapping)
df_reads_cell_good_bc["has_bc1"] = ~df_reads_cell_good_bc["bc1_correct"].isnull()
df_reads_cell_good_bc["bc2_collapse"] = df_reads_cell_good_bc.BC2.map(bc30_either_map)

# collapse UMIs
read_thresh = 4
cells_enough_reads = reads_per_cell.index[reads_per_cell >= read_thresh]
df_reads_enough = df_reads_cell_good_bc[
    df_reads_cell_good_bc.cell.isin(cells_enough_reads)
]
cell_list = [c[1] for c in df_reads_enough[["cell", "UB"]].groupby("cell")]
print(len(cell_list))

UMI_HAMMING = 2


def fix_umi(df_this_cell, umi_hamming=UMI_HAMMING):
    umi_vc = df_this_cell.UB.value_counts()
    umi_mapper = {}

    for i, umi in enumerate(umi_vc.index):
        # for each barcode go from n+1:end for
        if umi in umi_mapper:
            continue
        test_umis = umi_vc.index[i + 1 :].difference(umi_mapper.keys())
        for test in test_umis:
            if hamming(umi, test) <= umi_hamming:
                umi_mapper[test] = umi
    # add back original ones
    for umi in umi_vc.index:
        if umi not in umi_mapper:
            umi_mapper[umi] = umi
    return df_this_cell.UB.map(umi_mapper)


update_umi = parallel_process(cell_list, fix_umi, front_num=0)
df_reads_enough_umi = df_reads_enough.join(
    pd.concat(update_umi).rename("UB_collapse").sort_index()
)


# get consensus barcode for each UMI
cell_list = [c for c in df_reads_enough_umi.groupby("cell")]
bc_counts = parallel_process(cell_list, get_cell_barcodes, n_jobs=cores, front_num=0)
all_bc_counts = parallel_process(cell_list, cell_barcodes_frac, n_jobs=cores, front_num=0)

# summarize barcodes and umis for each cell
all_bc_flat = flatten(all_bc_counts)
df_cell_bc_counts = pd.DataFrame(
    all_bc_flat,
    columns=[
        "cell",
        "umi",
        "umi_rank",
        "BC1",
        "BC2",
        "bc_reads",
        "umi_reads",
        "n_above",
        "ratio_top",
        "frac_all",
        "frac_5",
        "frac_10",
    ],
)

df_cell_bc_counts["frac"] = (
    df_cell_bc_counts["bc_reads"] / df_cell_bc_counts["umi_reads"]
)
with_bc = [b for b in bc_counts if len(b) > 0]
print(len(with_bc))

key_list = []
bc_list = []
for v in with_bc:
    for k, vv in v.items():
        key_list.append(k)
        bc_list.append(vv)

df_with_bc = pd.DataFrame(
    bc_list, columns=["BC1", "BC2", "bc_reads", "umi_reads"]
).join(pd.DataFrame(key_list, columns=["cell", "umi"]))

df_all_bc = df_cell_bc_counts
df_bc_merge = df_all_bc.merge(
    df_with_bc, on=["cell", "umi", "BC2"], how="left", suffixes=(None, "_t")
)
name_base = MOUSE + "_"
df_bc_merge["cell_name"] = name_base + df_bc_merge["cell"]
assert df_bc_merge["cell_name"].isin(adata.obs_names).all()
CROUND = MOUSE.split("_")[1]

df_reads_cell_good_bc.reset_index(drop=True).to_parquet(
    save_data_folder / f"df_reads_cell_good_bc_{CROUND}.parquet"
)
df_all_bc.to_csv(save_data_folder / f"df_all_bc_umi_thresh_4_{CROUND}.csv")
df_with_bc.to_csv(save_data_folder / f"df_with_bc_umi_4_frac_075_{CROUND}.csv")
df_bc_merge.to_parquet(save_data_folder / f"df_bc_merge_umi_thresh_4_{CROUND}.parquet")