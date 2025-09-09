import os
import gzip
import pandas as pd
from polars import col
import scanpy as sc
from scipy.sparse import csr_matrix
from seaborn import heatmap
import squidpy as sq
import matplotlib.pyplot as plt
from pathlib import Path
ad = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/immune.scanvi.h5ad')
ad =ad.copy()
ad.obs['immune_celltype_coarse'].cat.categories
Macro = ad[ad.obs['immune_celltype_coarse']=='Macro']
Macro.obs['immune_celltype'].cat.categories

sc.tl.pca(Macro)

import scanpy.external as sce
sce.pp.harmony_integrate(Macro, 'Sample')
sc.pp.neighbors(Macro,use_rep='X_pca_harmony')
sc.tl.umap(Macro)
sc.tl.leiden(Macro,resolution=0.27)
sc.pl.umap(Macro,color=['leiden'],edges=False, legend_loc='on data')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/immune/Macro_leiden.pdf", dpi=1200, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/immune/Macro_leiden.png", dpi=1200, bbox_inches='tight')
sc.pl.umap(Macro,color='immune_celltype')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/immune/Macro_celltype.pdf", dpi=1200, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/immune/Macro_celltype.png", dpi=1200, bbox_inches='tight')
map = {
    '0':'0',
    '1':'1',
    '2':'2',
    '3':'3',
    '4':'4',
    '5':'5',
    '6':'6' 
}
sc.tl.rank_genes_groups(Macro, 'immune_celltype', method='wilcoxon')
sc.pl.rank_genes_groups(Macro, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/immune/Macro_rank_genes_groupsimmune_celltype.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/immune/Macro_rank_genes_groupsimmune_celltype.png",dpi = 300, bbox_inches='tight')
sc.pl.umap(Macro,color = marker_genes, use_raw=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/immune/Macro_genes.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/immune/Macro_genes.png",dpi = 300, bbox_inches='tight')

# 提取基因表达值（若存在）
import numpy as np
# 检查整个表达矩阵
if np.isnan(Macro.X).any():
    print("WARNING: 表达矩阵中存在 NaN 值")
    # 全局填充 NaN
    Macro.X = np.nan_to_num(Macro.X, nan=0.0)

marker_genes = [
                'CCL18',     #Artery·
                'APOE',     #Artery·
                'APOC1',
                'CHI3L1',     #Artery·
                'S100A6',     #Artery·
                'S100A10',
                'CXCL3',     #Artery·
                'CCL3',     #Artery·
                'CD83',
                'CXCL9',     #B
                'CXCL10',     #B
                'GBP1',
                'SELENOP',    #CD4
                'CD163',    #CD4
                'F13A1',
                'SPP1',
                'LDHA',
                'VCAN'
]
# 执行矩阵可视化
sc.pl.matrixplot(
    Macro, 
    marker_genes, 
    'immune_celltype', 
    standard_scale='var',
    colorbar_title='column scaled\nexpression'
)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/immune/Macro_heatmap.pdf",dpi = 300, bbox_inches='tight')
