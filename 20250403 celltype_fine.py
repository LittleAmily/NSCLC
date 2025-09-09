import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import anndata as ad
import scanpy as sc
import seaborn as sns
adata = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/un_integreted_data.h5ad')
adata = adata.copy()
adata
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
sns.jointplot(
    data=adata.obs,
    x="total_counts",
    y="n_genes_by_counts",
    kind="hex",
)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/beforeQC-1.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/beforeQC-1.png',dpi = 300, bbox_inches='tight')
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.1,
             groupby = 'Dataset', rotation= 90)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/beforeQC-2.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/beforeQC-2.png',dpi = 300, bbox_inches='tight')
sc.pl.violin(adata,['pct_counts_mt', 'pct_counts_hb', 'pct_counts_ribo'], jitter=0.1,
             groupby = 'Dataset', rotation= 90)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/beforeQC-3.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/beforeQC-3.png',dpi = 300, bbox_inches='tight')
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs['n_genes_by_counts'] <= 4000, :]
adata = adata[adata.obs['total_counts'] <= 20000, :]
adata = adata[adata.obs['pct_counts_mt'] <= 25, :]
adata = adata[adata.obs['pct_counts_hb'] <= 10, :]
sns.jointplot(
    data=adata.obs,
    x="total_counts",
    y="n_genes_by_counts",
    kind="hex",
)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/afterQC-1.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/afterQC-1.png',dpi = 300, bbox_inches='tight')
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.1,
             groupby = 'Dataset', rotation= 90)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/afterQC-2.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/afterQC-2.png',dpi = 300, bbox_inches='tight')
sc.pl.violin(adata,['pct_counts_mt', 'pct_counts_hb', 'pct_counts_ribo'], jitter=0.1,
             groupby = 'Dataset', rotation= 90)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/afterQC-3.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/afterQC-3.png',dpi = 300, bbox_inches='tight')
# 校正测序深度
sc.pp.normalize_total(adata, target_sum=1e4)
# 表达值log转化
sc.pp.log1p(adata)
# 备份完整矩阵
adata.raw = adata
# 计算高变基因
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
# 高变基因展示
sc.pl.highly_variable_genes(adata)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/highly_variable_genes.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/highly_variable_genes.png',dpi = 300, bbox_inches='tight')
# 剔除技术变异的影响
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
# scale data, clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/pca_variance_ratio.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/pca_variance_ratio.png',dpi = 300, bbox_inches='tight')
sc.pp.neighbors(adata, use_rep='X_scanvi_fix_linear')
sc.tl.umap(adata,n_components=20)
sc.pl.umap(adata, color=['C_scANVI'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_C_scANVI.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_C_scANVI.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=[ '_prediction'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_prediction.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_prediction.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Celltype'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Celltype.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Celltype.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Dataset'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Dataset.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Dataset.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Platform'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Platform.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Platform.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Drug'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Drug.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Drug.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Timepoint'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Timepoint.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Timepoint.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Tissue'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Tissue.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Tissue.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Origin'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Origin.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_Origin.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['batch'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_batch.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/SCANVI/umap_batch.png',dpi = 300, bbox_inches='tight')
# 恢复log1p的表达矩阵
adata = adata.raw.to_adata()
adata.write_h5ad("/home/data/sdzl14/NSCLC/zong/adata_log.h5ad", compression="gzip")
adata = sc.read_h5ad("/home/data/sdzl14/NSCLC/zong/adata_log.h5ad")
adata = adata.copy()
adata.X.max()
raw = sc.read_h5ad("/home/data/sdzl14/NSCLC/zong.first_inter.h5ad")
raw.X.max()
# 对齐 raw 和 adata
subset = raw[raw.obs_names.isin(adata.obs_names), raw.var_names.isin(adata.var_names)]
# 将 raw 的 counts 数据添加到 adata.layers
adata.layers['counts'] = subset.X
adata.write_h5ad("/home/data/sdzl14/NSCLC/zong/adata_log.h5ad", compression="gzip")
sc.pp.neighbors(adata, use_rep='X_scanvi_fix',method='umap')
sc.tl.umap(adata)
sc.pl.umap(adata, color=['C_scANVI'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_C_scANVI.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_C_scANVI.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['C_scANVI'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_C_scANVI.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_C_scANVI.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=[ '_prediction'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_prediction.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_prediction.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Celltype'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Celltype.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Celltype.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Dataset'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Dataset.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Dataset.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Platform'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Platform.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Platform.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Drug'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Drug.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Drug.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Timepoint'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Timepoint.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Timepoint.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Tissue'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Tissue.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Tissue.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['Origin'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Origin.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_Origin.png',dpi = 300, bbox_inches='tight')
sc.pl.umap(adata, color=['batch'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_batch.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_batch.png',dpi = 300, bbox_inches='tight')
adata = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/third.scANVI_adata.h5ad')