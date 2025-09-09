import numpy as np

import cellrank as cr
import scanpy as sc
import scvelo as scv
adata = sc.read_h5ad("/home/data/sdzl14/NSCLC/zong/epi.h5ad")
adata = adata.copy()
scv.pp.filter_and_normalize(
    adata, min_shared_counts=20, n_top_genes=2000, subset_highly_variable=False
)

sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)