import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import anndata as ad
import scanpy as sc
import seaborn as sns
import scanorama
adata=sc.read_h5ad("/home/data/sdzl14/NSCLC/zong/adata_log.h5ad")
adata = adata.copy()
sc.pp.neighbors(adata, use_rep='X_scanvi_fix',method='umap')
sc.tl.umap(adata)
sc.pl.umap(adata, color=['C_scANVI'], edges=False)
sc.pl.umap(adata, color=['Celltype'], edges=False)
adata.var['highly_variable'].value_counts()

macrophage = adata[adata.obs['C_scANVI'] == 'Macrophage']
sc.pl.umap(adata, color=['Dataset'], edges=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_C_scANVI_Dataset.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_C_scANVI_Dataset.png',dpi = 300, bbox_inches='tight')
adata.obs['_prediction'].value_counts()
immune = ['CD4+T','CD8+T','Macrophage','NK','Neutrophil','Memory B cell','Plasma','Monocyte','Treg','DC','Naive B cell','Mast','pDC','GC B cell','ILC']
immune = adata[adata.obs['_prediction'].isin(immune)]
sc.pl.umap(immune,color = '_prediction',)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/umap_C_scANVI_immune.pdf',dpi = 300, bbox_inches='tight')
plt.savefig ('/home/data/sdzl14/NSCLC/zong/fig/umap_C_scANVI_immune.png',dpi = 300, bbox_inches='tight')
ad5= sc.read_h5ad('/home/data/sdzl14/NSCLC/Ansuman_Satpathy_2023/Ansuman_Satpathy_2023_filtered.h5ad')
ad5
ad1 = sc.read_h5ad('/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_annotated.h5ad')
ad1 = ad1.copy()
ad2 = sc.read_h5ad('/home/data/sdzl14/NSCLC/Peng_Zhang_2024/scRNA_annotated.h5ad')
ad2 = ad2.copy()
ad3 = sc.read_h5ad('/home/data/sdzl14/NSCLC/SDB_integrated/SDB_integrated.h5ad')
ad3 = ad3.copy()

ad4 = sc.read_h5ad('/home/data/sdzl14/NSCLC/Tagore_S_2025/Tagore_S_2025_integrated_data.h5ad')
ad4  = ad4.copy()
ad1.obs['celltype_fine'].cat.categories
immune1 = ['CD4+T', 'CD8+T','DC','Macrophage',
       'Mast cells', 'Memory B cells', 'Monocytes', 'NK cells',
       'Naive B cells', 'Neutrophils', 'Plasma cells', 'pDC']
immune1 =ad1[ad1.obs['celltype_fine'].isin(immune1)] 
ad2.obs['minor'].cat.categories
immune2 = ['CCL3L1+ cDC2',
       'CD1A+CD207+ cDC2', 'CD1C+ mregDC', 'CD4_CXCL13', 'CD4_NR4A2',
       'CD4_TCF7', 'CD8_GZMK', 'CD8_HAVCR2', 'CD8_IFNG','Cycling B', 'Cycling T','GC B', 'ILC','MB_FCRL4', 'MB_NR4A1',
       'MB_TXNIP', 'LTB+ cDC2','MRC1+ cDC2', 'MRC1+IL1B+ cDC2','Macro_CCL18', 'Macro_CHI3L1', 'Macro_CXCL3', 'Macro_CXCL9',
       'Macro_SELENOP', 'Macro_SPP1', 'AM','Mast', 'Mono_FCGR3A', 'Mono_VEGFA',
       'NK_FCGR3A', 'NK_GNLY', 'Naive B', 'Neutrophil', 'Pericyte', 'Plasma','Treg','cDC1', 'mregDC','pDC']
immune2 =ad2[ad2.obs['minor'].isin(immune2)]
ad3.obs['cell_type'].cat.categories
immune3 = ['B cell',
       'B cell dividing','DC mature','Macrophage', 'Macrophage alveolar',
       'Mast cell', 'Monocyte classical',
       'Monocyte non-classical', 'NK cell', 'NK cell dividing', 'Neutrophils',
       'Pericyte', 'Plasma cell', 'Plasma cell dividing', 'T cell CD4',
       'T cell CD4 dividing', 'T cell CD8 activated', 'T cell CD8 dividing',
       'T cell CD8 effector memory', 'T cell CD8 naive',
       'T cell CD8 terminally exhausted', 'T cell NK-like',
       'T cell regulatory','cDC1', 'cDC2', 'myeloid dividing',
       'pDC']
immune3 = ad3[ad3.obs['cell_type'].isin(immune3)]
ad4.obs['tumor_nontumor_finer'].cat.categories
immune4 = ['B_cells','MDMs/Monocytes/Macrophage 1 ', 'Microglia', 'Myeloid (Others)',
       'Neutrophils','T_cells','Alveolar Macrophages'
]
immune4 = ad4[ad4.obs['tumor_nontumor_finer'].isin(immune4)]
immune5 = ad5
immune = immune1.concatenate(immune2,immune3,immune4,immune5,join="inner")
keep_obs_columns = ['Sample', 'Patient', 'Celltype', 'Dataset', 'Platform', 
                    'Pathtype', 'Drug', 'Timepoint', 'Tissue', 'Origin','tumor_nontumor_finer','minor','cell_type','celltype_fine']
immune.obs = immune.obs[keep_obs_columns].copy()
# 清理 var 的元数据（保留基因名，删除所有列）
immune.var = pd.DataFrame(index=immune.var.index)  # 仅保留基因名，删除所有列

# 删除 obsm 中的嵌入数据
immune.obsm = {}

# 删除 layers 中的额外数据层
immune.layers = {}
immune.uns = {}
immune.obsp = {}
immune.write_h5ad("/home/data/sdzl14/NSCLC/zong/immune.h5ad")
immune_scanvi = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/immune.scanvi.h5ad')
immune_scanvi = immune_scanvi.copy()
immune = sc.read_h5ad("/home/data/sdzl14/NSCLC/zong/immune.h5ad")
immune = immune.copy()
print(immune_scanvi.obs['Dataset'].value_counts())
print(immune.obs['Dataset'].value_counts())
# 将多个数据集名称放入列表中
selected_datasets = ['Ansuman_Satpathy_2023', 'Goveia_Carmeliet_2020', 
       'He_Fan_2021', 'Kim_Lee_2020', 'Lambrechts_Thienpont_2018_6149v1',
       'Lambrechts_Thienpont_2018_6149v2', 'Lambrechts_Thienpont_2018_6653',
       'Laughney_Massague_2020', 'Maynard_Bivona_2020', 'Peng_Zhang_2024',
       'Tagore_S_2025', 'Travaglini_Krasnow_2020', 'UKIM-V',
       'Vieira_Teichmann_2019']
immune = immune[immune.obs['Dataset'].isin(selected_datasets)]

print(immune.obs['Dataset'].value_counts())
immune.layers['counts'] = immune.X.copy()

immune.obsm['X_scanvi_fix'] = immune_scanvi.obsm['X_scanvi_fix'].copy()
immune.obsm['X_scanvi_fix_linear'] = immune_scanvi.obsm['X_scanvi_fix_linear'].copy()
immune.obsm['X_scanvi_no_fix'] = immune_scanvi.obsm['X_scanvi_no_fix'].copy()
immune.obsm['X_scvi'] = immune_scanvi.obsm['X_scvi'].copy()
immune.obsm['X_umap'] = immune_scanvi.obsm['X_umap'].copy()
immune.obsm['_latent'] = immune_scanvi.obsm['_latent'].copy()
print(immune_scanvi.obs['_prediction'].value_counts())
print(immune_scanvi.obs['_scvi_batch'].value_counts())
immune.obs['immune_celltype'] = immune_scanvi.obs['_prediction'].copy()
sc.pp.normalize_total(immune, target_sum=1e4)
sc.pp.log1p(immune)
sc.pp.highly_variable_genes(
    immune,
    flavor="seurat",
    n_top_genes=2000,

    inplace=True,batch_key="Dataset")

sc.pp.scale(immune, max_value=10)
sc.pp.neighbors(immune, use_rep='X_scanvi_fix_linear',n_neighbors=35)
sc.tl.umap(immune,min_dist=0.5)
sc.pl.umap(immune,color=['immune_celltype'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/immune_celltype.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/immune_celltype.pdf',dpi = 300, bbox_inches='tight')
sc.tl.tsne(immune,n_jobs=20,use_rep='X_scanvi_fix_linear')
sc.pl.tsne(immune,color=['immune_celltype','Dataset'],legend_loc='on data')
print(immune_scanvi.obs['Tissue'].value_counts())
sc.pl.umap(immune,color=['Tissue'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/Tissue.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/Tissue.pdf',dpi = 300, bbox_inches='tight')
sc.pl.umap(immune,color=['Origin'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/Origin.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/Origin.pdf',dpi = 300, bbox_inches='tight')
sc.pl.umap(immune,color=['Dataset'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/Dataset.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/Dataset.pdf',dpi = 300, bbox_inches='tight')
sc.tl.rank_genes_groups(immune, 'immune_celltype', method='wilcoxon')
sc.pl.rank_genes_groups(immune, n_genes=20, sharey=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/rank_genes_groups.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/rank_genes_groups.pdf',dpi = 300, bbox_inches='tight')