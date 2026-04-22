import joblib
import numpy as np
import pandas as pd
from numba import njit
from pynndescent import NNDescent
from scipy import sparse, stats
from tqdm import trange


def mmnorm(x):
    xi = x.min()
    xa = x.max()
    return (x - xi) / (xa - xi)


k_neigh = 20
beta = 0.1
beta_i = 1 - beta
n_iter = 1500

adata = ""

# KNN graph from scVI latent space
X_pca = adata.obsm["X_INP5_scVI"]
index = NNDescent(X_pca, metric="cosine", verbose=True, n_jobs=-1)
indices, distances = index.query(X_pca, k=k_neigh + 1)

# smoother matrix for KNN-smoothing
n_cells = X_pca.shape[0]
row = np.repeat(np.arange(n_cells), k_neigh)
col = indices[:, 1:].flatten()
data = np.ones_like(row) / k_neigh
a_mat = sparse.csr_matrix((data, (row, col)), shape=(n_cells, n_cells))
assert np.allclose(np.ones(n_cells), a_mat.sum(1).A.flatten())
assert ((a_mat > 0).sum(1) == k_neigh).all()
ind = sparse.eye(n_cells)
# weight by itself vs neighbors
smoother = beta * ind + beta_i * a_mat

# start from a single cell and run diffusion until converges
root_idx = np.where(
    (adata.obsm["celltype_prob"]["GBC"] == 1) & (adata.obs.source == "homecage")
)[0][-1]

pseudo_vals = np.zeros(n_cells)
pseudo_vals[root_idx] = 1
out = pseudo_vals.copy()
all_outs = np.zeros((n_cells, n_iter))

for i in trange(n_iter):
    out = mmnorm(smoother @ out)
    all_outs[:, i] = out

ranks = stats.rankdata(all_outs, axis=0)
rank_frac = 1 - ranks[:, -1] / ranks[:, -1].max()
adata.obs["pseudotime_rank"] = rank_frac
adata.obs["pseudotime_dist"] = all_outs[:, -1]

out_dict = {
    "root_idx": root_idx,
    "X_pca": X_pca,
    "index": index,
    "indices": indices,
    "obs": adata.obs,
    "k_neigh": k_neigh,
}
joblib.dump(out_dict, "INP5_pseudotime_smoother_index.p")


### example gene smoothing
# list of genes
to_smooth = []
adata_gene_subset = adata[:, to_smooth].copy()

orig = adata_gene_subset.X.A
out = orig.copy()
out[np.isnan(out)] = 0
# pick number of times to smooth
for i in trange(5):
    out = smoother @ out

adata_gene_subset.X = out.copy()
df_gene_smooth = pd.DataFrame(
    adata_gene_subset.X,
    index=adata_gene_subset.obs_names,
    columns=adata_gene_subset.var_names,
)
