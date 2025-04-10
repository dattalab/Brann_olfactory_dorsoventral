import gzip
import logging
import os
from pathlib import Path

import pandas as pd


def unique_path(directory, name_pattern):
    """Iterate over file names until find one that doesn't exists
    Args:
        directory (pathlib.Path)
        name_pattern (str): with formatter like "counts_no_double_umi_{:03d}.tsv.gz"
    """
    counter = 0
    while True:
        counter += 1
        path = directory / name_pattern.format(counter)
        if not path.exists():
            return path


def load_barcodes(bcfile):
    """Read 10x-genomices barcodes.tsv.gz file"""
    if ".gz" in bcfile.suffixes:
        f = gzip.open(bcfile, "rt")
    else:
        f = open(bcfile, "r")
    valid_bcset = set([s.rstrip("\n") for s in f.readlines()])
    return valid_bcset


def find_and_load_barcodes(bam_fn, bcfile=None):
    """walk down 10x outs dir and find bcfile"""
    if bcfile is not None:
        for dir in [x for x in bam_fn.parent.glob("filtered_*_bc_matr*") if x.is_dir()]:
            # check if bcfile in folder
            files = [x for x in dir.rglob("barcodes.tsv*") if x.is_file()]
            if len(files) > 0:
                bcfile = files[0]  # take first batch
                logging.debug(f"Found bcfile {bcfile} in {dir}")
                break
    if not bcfile or not bcfile.is_file():
        raise OSError(f"Could not find bcfile")
    else:
        return load_barcodes(bcfile)


def get_cores(verbose=False):
    """Get number of cores, checking for limits if on SLURM"""
    try:
        cores = int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
    except KeyError as ke:
        cores = os.cpu_count()
    if verbose:
        print(f"Found {cores} cpus")
    return cores


def filter_cell_reads(tags, current_cell, out_fn=None, exclude_multi_gene=True):
    """Removed double-counted UMIs from a list of reads.
    Follows 10x logic to ignore reads that mapped to a region with multiple genes (with semicolon)
    Removes umis where different reads each uniquely mapped to separate genes, then
    counts the number of UMIs for each gene and saves it in tidy data format (ala a sparse matrix)

    Args:
        tags (list of tuples): GN (gene) and UB (UMI) tag for each read
        current_cell (str): CB (cell barcode) for reads
        out_fn (pathlib.Path): tsv file to save UMI counts for each gene

    Returns:
        current_cell (cell): name of cell processed
    """
    df = pd.DataFrame(tags, columns=["GN", "UB"])
    if exclude_multi_gene:
        df_keep = df[~df.GN.str.contains(";", na=False)]
    else:
        df_keep = df
    # number of genes each UMI mapped to
    umi_counts = df_keep.groupby("UB").GN.nunique()
    # exclude UMIs that uniquely mapped to multiple genes
    df_to_count = df_keep[~df_keep.UB.isin(umi_counts.index[umi_counts > 1])]
    df_to_write = df_to_count.groupby("GN", as_index=False).UB.nunique()
    df_to_write["CB"] = current_cell.split("-")[0]
    df_to_write.to_csv(out_fn, sep="\t", mode="a", header=None, index=False)
    return current_cell