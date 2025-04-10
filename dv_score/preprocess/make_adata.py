import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.sparse import csr_matrix


def run(fn, outfn=None, prefix=None):
    out_folder = fn.parents[0]
    if prefix is None:
        prefix = fn.parents[1].name
    if not prefix.endswith("_"):
        prefix += "_"
    if outfn is None:
        outfn = Path(out_folder, prefix + fn.with_suffix("").stem).with_suffix(".h5ad")
    print(f"loading file {fn}")
    print(f"using prefix {prefix}")
    df = pd.read_csv(fn, sep="\t", header=None, index_col=None)
    if "cell" in df.columns:
        df = df.rename({"count": "n"}, axis="columns")
    else:
        df.columns = ["gene", "n", "cell"]
    print(f"Loaded, shape = {df.shape}")
    for col in ("gene", "cell"):
        print(col)
        df[col] = df[col].astype("category")
    # rows are cells
    # cols are genes
    mat = csr_matrix(
        (df.n.values, (df.cell.cat.codes, df.gene.cat.codes)), dtype="float32"
    )
    adata = AnnData(mat)
    adata.obs_names = prefix + df.cell.cat.categories
    adata.var_names = df.gene.cat.categories
    print(adata)
    adata.var["mito"] = adata.var_names.str.contains("^mt-")
    adata.var["ribo"] = adata.var_names.str.contains("^Rp[ls]")
    adata.obs["orig_ident"] = prefix.rstrip("_")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mito", "ribo"], percent_top=None, inplace=True
    )
    adata.obs["n_counts"] = adata.obs["total_counts"]
    adata.obs["n_genes"] = adata.obs["n_genes_by_counts"]
    print(adata.obs.n_counts.describe())
    print(f"Saving to outfn {outfn}")
    adata.write_h5ad(Path(outfn))
    print("Saved.")


def parse_args():
    # parse args
    parser = argparse.ArgumentParser()
    path_type = lambda p: Path(p).absolute()
    parser.add_argument(
        "-f",
        "--fn",
        dest="fn",
        help="filename",
        type=path_type,
        default="counts_no_double_umi_001.tsv",
    )
    parser.add_argument(
        "-o", "--out", dest="out", help="output filename", type=str, default=None
    )
    parser.add_argument(
        "-p",
        "--pre",
        dest="pre",
        help="prefix for cell barcodes",
        type=str,
        default=None,
    )

    return parser.parse_args()

def main():
    args = parse_args()
    run(args.fn, args.out, args.pre)


if __name__ == "__main__":
    main()
