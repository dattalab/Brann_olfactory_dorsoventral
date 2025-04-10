import argparse
import logging
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from pathlib import Path

import pysam

from dv_score.preprocess import util

logging.basicConfig(
    stream=sys.stdout,
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO,
)


def parse_args():
    # parse args
    parser = argparse.ArgumentParser()
    path_type = lambda p: Path(p).absolute()
    parser.add_argument(
        "-i",
        "--in_bam",
        dest="in_bam",
        help="input cell_sorted bam filename. \
            Output of `samtools view -b -q 255 -e '(![MM] && [GN]) && [UB]' -D CB:barcodes.tsv` possorted_genome_bam.bam \
            and `samtools sort -t CB possorted_genome_bam.unique.GN.barcode.bam`",
        type=path_type,
        default="cell_sorted_possorted_genome_bam.unique.bam",
    )
    parser.add_argument(
        "-b",
        "--bcfile",
        dest="bcfile",
        help="barcodes tsv file",
        type=path_type,
        default=None,
    )
    parser.add_argument(
        "-o",
        "--out_fn",
        dest="out_fn",
        help="output filename",
        type=path_type,
        default=None,
    )
    parser.add_argument(
        "-c",
        "--cores",
        dest="cores",
        help="number of processors",
        type=int,
        default=15,
    )
    return parser.parse_args()


def run(bam_fn, bcfile, out_fn, n_cores, print_every=5_000_000):
    if not bam_fn.exists():
        raise ValueError(f"bam file {bam_fn} does not exists.")
    if out_fn is None:
        out_fn = util.unique_path(bam_fn.parent, "counts_no_double_umi_{:03d}.tsv")
    elif out_fn.exists():
        raise ValueError(f"out_fn {out_fn} already exists.")
    logging.info(f"Saving file to {out_fn}")
    valid_bcset = util.find_and_load_barcodes(bam_fn, bcfile)
    total_cell_number = len(valid_bcset)
    cores = min(n_cores, util.get_cores())
    logging.info(f"Using {cores} cores")

    current_cell = None
    tags = []
    n_cell = 0
    futures = []
    process = partial(util.filter_cell_reads, out_fn=out_fn)

    st = time.time()
    with pysam.AlignmentFile(bam_fn) as fin:
        with ProcessPoolExecutor(max_workers=cores) as executor:
            for i, read in enumerate(fin):
                if i % print_every == 0:
                    et = (time.time() - st) / 60  # in minutes
                    rt = (et * total_cell_number / (n_cell + 1)) - et
                    logging.info(
                        f"Filtered {n_cell:05d} of {total_cell_number} cells ({n_cell / total_cell_number * 100 :06.3f} %)"
                        + f" in first {i // 1000000} million reads in {et:04.2f} minutes...{rt:04.2f} minutes left."
                    )
                this_cell = read.get_tag("CB")
                if this_cell != current_cell:
                    if n_cell > 0:
                        future = executor.submit(process, tags, current_cell)
                        futures.append(future)
                    n_cell += 1
                    current_cell = this_cell
                    tags = []
                tags.append((read.get_tag("GN"), read.get_tag("UB")))

            # finish reading and process last cell
            logging.info(f"Processing cell {n_cell} with CB: {current_cell} and {len(tags)} tags")
            future = executor.submit(process, tags, current_cell)
            futures.append(future)

    logging.info("Finished reading bam file, checking all cells were saved")
    cells_done = []
    for future in as_completed(futures):
        try:
            cells_done.append(future.result())
        except Exception as e:
            cells_done.append(e)
    cells_done = set(cells_done)
    if cells_done == valid_bcset:
        logging.info("Found all cells.")
        logging.info("Removing writing permissions to the ouput tsv file")
        logging.info(oct(os.stat(out_fn).st_mode))
        # prevent future appending
        os.chmod(out_fn, 0o444)
        logging.info(oct(os.stat(out_fn).st_mode))
    else:
        raise ValueError("Did not find all cells, erroring.")


def main():
    args = parse_args()
    run(args.in_bam, args.bcfile, args.out_fn, args.cores)


if __name__ == "__main__":
    main()
