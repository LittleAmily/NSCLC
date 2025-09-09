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
ad1 = sc.read_h5ad('/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_annotated.h5ad')
ad1 = ad1.copy()
ad2 = sc.read_h5ad('/home/data/sdzl14/NSCLC/Peng_Zhang_2024/scRNA_annotated.h5ad')
ad2 = ad2.copy()
ad2.obs['PathType'].value_counts()
old_to_new =  {
'AT1':'Normal_Epi', 
'AT2':'Normal_Epi', 
'CD4+T':'CD4+T',
'CD8+T':'CD8+T', 
'Cancer_P1/2/3':'Malignant', 
'Cancer_P4':'Malignant', 
'Cancer_P5':'Malignant', 
'Ciliated cells':'Normal_Epi', 
'Club cells':'Normal_Epi', 
'Cycling cells':'Cycling cell', 
'DC':'DC',
'Endothelial cells':'Endothelium', 
'Epithelial cells':'Normal_Epi', 
'Fibroblasts':'Fibroblast', 
'Macrophage':'Macrophage',
'Mast cells':'Mast', 
'Memory B cells':'Memory B cell', 
'Monocytes':'Monocyte', 
'NK cells':'NK',
'Naive B cells':'Naive B cell', 
'Neutrophils':'Neutrophil', 
'Plasma cells':'Plasma', 
'pDC':'pDC'}
ad1.obs['Celltype_fine'] = (
    ad1.obs['celltype_fine']
    .map(old_to_new)
    .astype('category')
)
old_to_new ={
    'Artery':'Endothelium', 
    'CAF':'Fibroblast', 
    'CD4 T':'CD4+T', 
    'CD8 T':'CD8+T', 
    'Capillary':'Endothelium', 
    'Cycling B':'Cycling cell',
    'Cycling T':'Cycling cell', 
    'GC B':'GC B cell', 
    'ILC':'ILC', 
    'Lymphatic EC':'Endothelium', 
    'Macrophage':'Macrophage', 
    'Malignant':'Malignant',
    'Mast':'Mast', 
    'Memory B':'Memory B cell', 
    'Monocyte':'Monocyte', 
    'NK':'NK', 
    'Naive B':'Naive B cell',
    'Neutrophil':'Neutrophil',
    'Normal.Epi':'Normal_Epi', 
    'Pericyte':'Endothelium', 
    'Plasma':'Plasma', 
    'SMC':'Endothelium', 
    'Tip':'Endothelium', 
    'Treg':'Treg', 
    'Venule':'Endothelium',
    'cDC1':'DC', 
    'cDC2':'DC', 
    'mregDC':'DC', 
    'pDC':'pDC'
}
ad2.obs['Celltype_fine'] = (
    ad2.obs['major']
    .map(old_to_new)
    .astype('category')
)
ad3 = sc.read_h5ad('/home/data/sdzl14/NSCLC/SDB_integrated/SDB_integrated.h5ad')
ad3 = ad3.copy()

ad4 = sc.read_h5ad('/home/data/sdzl14/NSCLC/Tagore_S_2025/Tagore_S_2025_integrated_data.h5ad')
ad4  = ad4.copy()
keep_obs_columns = ['Sample', 'Patient', 'Celltype', 'Dataset', 'Platform', 
                    'Pathtype', 'Drug', 'Timepoint', 'Tissue', 'Origin','Celltype_fine']
ad2.obs['Pathtype'] = ad2.obs['PathType']
ad2.obs = ad2.obs[keep_obs_columns].copy()

# 清理 var 的元数据（保留基因名，删除所有列）
ad2.var = pd.DataFrame(index=ad2.var.index)  # 仅保留基因名，删除所有列

# 删除 obsm 中的嵌入数据
ad2.obsm = {}

# 删除 layers 中的额外数据层
ad2.layers = {}
ad2.uns = {}
ad2.obsp = {}
ad2

ad1.obs['Dataset'] = 'Peilin_Wang_2025'
ad1.obs = ad1.obs[keep_obs_columns].copy()

# 清理 var 的元数据（保留基因名，删除所有列）
ad1.var = pd.DataFrame(index=ad1.var.index)  # 仅保留基因名，删除所有列

# 删除 obsm 中的嵌入数据
ad1.obsm = {}

# 删除 layers 中的额外数据层
ad1.layers = {}
ad1.uns = {}
ad1.obsp = {}
ad1

keep_obs_columns = ['Sample', 'Patient', 'Celltype', 'Dataset', 'Platform', 
                    'Pathtype', 'Drug', 'Timepoint', 'Tissue', 'Origin','cell_type_major']

ad3.obs = ad3.obs[keep_obs_columns].copy()

# 清理 var 的元数据（保留基因名，删除所有列）
ad3.var = pd.DataFrame(index=ad3.var.index)  # 仅保留基因名，删除所有列

# 删除 obsm 中的嵌入数据
ad3.obsm = {}

# 删除 layers 中的额外数据层
ad3.layers = {}
ad3.uns = {}
ad3.obsp = {}
ad3

keep_obs_columns = ['Sample', 'Patient', 'Celltype', 'Dataset', 'Platform', 
                    'Pathtype', 'Drug', 'Timepoint', 'Tissue', 'Origin','tumor_nontumor_finer']
ad4.obs = ad4.obs[keep_obs_columns].copy()

# 清理 var 的元数据（保留基因名，删除所有列）
ad4.var = pd.DataFrame(index=ad4.var.index)  # 仅保留基因名，删除所有列

# 删除 obsm 中的嵌入数据
ad4.obsm = {}

# 删除 layers 中的额外数据层
ad4.layers = {}
ad4.uns = {}
ad4.obsp = {}
ad4
new_var_names = pd.read_csv("/home/data/sdzl14/NSCLC/Peng_Zhang_2024/rownames.txt", usecols=[0]).iloc[:,0].values

new_var_names
ad2.var_names = new_var_names
ad1.obs['celltype_main'].cat.categories
epi1 = ['AT1', 'AT2','Cancer_P1/2/3', 'Cancer_P4',
       'Cancer_P5','Ciliated cells', 'Club cells','Epithelial cells']
old_to_new = {
    'AT1':'AT1', 
    'AT2':'AT2',
    'Cancer_P1/2/3':'Tumor', 
    'Cancer_P4':'Tumor',
    'Cancer_P5':'Tumor',
    'Ciliated cells':'Ciliated cell', 
    'Club cells':'Club cell',
    'Epithelial cells':'Epithelial cell'}
epi1 =ad1[ad1.obs['celltype_fine'].isin(epi1)]
epi1.obs['celltype_fine'] = (
    epi1.obs['celltype_fine']
    .map(old_to_new)
    .astype('category')
)
ad2.obs['minor'].cat.categories
epi2 = ['E0_AT2',
       'E1_Malig', 'E2_Malig', 'E3_Malig', 'E4_Malig', 'E5_Malig', 'E6_Mucous',
       'E7_Malig', 'E8_Malig', 'E8_Normal', 'E9_Malig', 'E10_Malig',
       'E11_Malig', 'E12_Malig', 'E13_Malig', 'E14_Basal', 'E15_Ciliated',
       'E16_Club', 'E17_Serous', 'E18_Malig', 'E19_Malig', 'E19_Normal']
old_to_new = {
'E0_AT2':'AT2',
       'E1_Malig':'Tumor', 'E2_Malig':'Tumor', 'E3_Malig':'Tumor', 'E4_Malig':'Tumor', 'E5_Malig':'Tumor', 'E6_Mucous':'Tumor',
       'E7_Malig':'Tumor', 'E8_Malig':'Tumor', 'E8_Normal':'Tumor', 'E9_Malig':'Tumor', 'E10_Malig':'Tumor',
       'E11_Malig':'Tumor', 'E12_Malig':'Tumor', 'E13_Malig':'Tumor', 'E14_Basal':'Basal cell', 'E15_Ciliated':'Ciliated cell',
       'E16_Club':'Club cell', 'E17_Serous':'Serous cell', 'E18_Malig':'Tumor', 'E19_Malig':'Tumor', 'E19_Normal':'Epithelial cell'}
epi2 = ad2[ad2.obs['minor'].isin(epi2)]
epi2.obs['celltype_fine'] = (
    epi2.obs['minor']
    .map(old_to_new)
    .astype('category')
)
ad3.obs['cell_type'].cat.categories
epi3 = ['Alveolar cell type 1', 'Alveolar cell type 2',  'Ciliated', 'Club',
       'ROS1+ healthy epithelial', 'Tumor cells', 'transitional club/AT2']
old_to_new = {
'Alveolar cell type 1':'AT1', 'Alveolar cell type 2':'AT2',  'Ciliated':'Ciliated cell', 'Club':'Club cell',
       'ROS1+ healthy epithelial':'Epithelial cell', 'Tumor cells':'Tumor', 'transitional club/AT2':'AT2'}
epi3 = ad3[ad3.obs['cell_type'].isin(epi3)]
epi3.obs['celltype_fine'] = (
    epi3.obs['cell_type']
    .map(old_to_new)
    .astype('category')
)
ad4.obs['tumor_nontumor_finer'].cat.categories
epi4 = [ 'Tumor', 'Tumor (Not Tumor)','Other']
other = ad4[ad4.obs['tumor_nontumor_finer']== 'Other']
epi4 = ad4[ad4.obs['tumor_nontumor_finer'].isin(epi4)]
epi = epi1.concatenate(epi2, epi3, epi4, join= 'inner')
keep_obs_columns = ['Sample', 'Patient', 'Celltype', 'Dataset', 'Platform', 
                    'Pathtype', 'Drug', 'Timepoint', 'Tissue', 'Origin','tumor_nontumor_finer','minor','cell_type','celltype_fine']
epi.obs = epi.obs[keep_obs_columns].copy()
# 清理 var 的元数据（保留基因名，删除所有列）
epi.var = pd.DataFrame(index=epi.var.index)  # 仅保留基因名，删除所有列

# 删除 obsm 中的嵌入数据
epi.obsm = {}

# 删除 layers 中的额外数据层
epi.layers = {}
epi.uns = {}
epi.obsp = {}
epi.write_h5ad("/home/data/sdzl14/NSCLC/zong/epi.h5ad")
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

epi_scanvi = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/epi.scanvi.h5ad')
epi_scanvi = epi_scanvi.copy()
epi.obs_names
epi = sc.read_h5ad("/home/data/sdzl14/NSCLC/zong/epi.h5ad")
epi = epi.copy()
ad2.obs_names
print(epi_scanvi.obs['Dataset'].value_counts())
print(epi.obs['Dataset'].value_counts())
# 将多个数据集名称放入列表中
selected_datasets = ['Goveia_Carmeliet_2020', 'He_Fan_2021', 'Kim_Lee_2020',
       'Lambrechts_Thienpont_2018_6149v1', 'Lambrechts_Thienpont_2018_6149v2',
       'Lambrechts_Thienpont_2018_6653', 'Laughney_Massague_2020',
       'Maynard_Bivona_2020', 'Peng_Zhang_2024', 'Tagore_S_2025',
       'Travaglini_Krasnow_2020', 'Vieira_Teichmann_2019']
epi = epi[epi.obs['Dataset'].isin(selected_datasets)]

print(epi.obs['Dataset'].value_counts())
epi.layers['counts'] = epi.X.copy()

epi.obsm['X_scanvi_fix'] = epi_scanvi.obsm['X_scanvi_fix'].copy()
epi.obsm['X_scanvi_fix_linear'] = epi_scanvi.obsm['X_scanvi_fix_linear'].copy()
epi.obsm['X_scanvi_no_fix'] = epi_scanvi.obsm['X_scanvi_no_fix'].copy()
epi.obsm['X_scvi'] = epi_scanvi.obsm['X_scvi'].copy()
epi.obsm['X_umap'] = epi_scanvi.obsm['X_umap'].copy()
epi.obsm['_latent'] = epi_scanvi.obsm['_latent'].copy()
print(epi_scanvi.obs['_prediction'].value_counts())
print(epi_scanvi.obs['_scvi_batch'].value_counts())
epi.obs['Epithelium_type'] = epi_scanvi.obs['_prediction'].copy()
epi.obs['Pathtype'].value_counts()
sc.pp.normalize_total(epi, target_sum=1e4)
sc.pp.log1p(epi)
sc.pp.highly_variable_genes(
    epi,
    flavor="seurat",
    n_top_genes=2000,

    inplace=True,batch_key="Dataset")

sc.pp.scale(epi, max_value=10)
sc.pp.neighbors(epi, use_rep='X_scanvi_fix_linear',n_neighbors=35)
sc.tl.umap(epi,min_dist=0.5)
sc.pl.umap(epi,color=['Epithelium_type'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Epithelium_type.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Epithelium_type.pdf',dpi = 300, bbox_inches='tight')
sc.tl.tsne(epi,n_jobs=20,use_rep='X_scanvi_fix_linear')
sc.pl.tsne(epi,color=['Epithelium_type','Dataset'],legend_loc='on data')
print(epi_scanvi.obs['Tissue'].value_counts())
sc.pl.umap(epi,color=['Tissue'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Tissue.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Tissue.pdf',dpi = 300, bbox_inches='tight')
sc.pl.umap(epi,color=['Origin'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Origin.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Origin.pdf',dpi = 300, bbox_inches='tight')
sc.pl.umap(epi,color=['Dataset'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Dataset.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Dataset.pdf',dpi = 300, bbox_inches='tight')
sc.tl.rank_genes_groups(epi, 'Epithelium_type', method='wilcoxon')
sc.pl.rank_genes_groups(epi, n_genes=20, sharey=False)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/rank_genes_groups.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/rank_genes_groups.pdf',dpi = 300, bbox_inches='tight')
marker_genes = ['CAV1',   #AT1
                'SFTPA1',     #AT2
                'KRT17',    #Basal
                'TPPP3',   #Ciliated
                'SDK1',     #Club
                'CNTN5',  #Normal
                'SYNE2',   #Serous
                'GAPDH'    #Tumor

]
sc.pl.umap(epi, color=marker_genes)
plt.savefig('/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/20250107sdb/13_umapleiden_marker_genes.pdf',dpi = 300, bbox_inches='tight')
sc.pl.matrixplot(epi, marker_genes, 'Epithelium_type', standard_scale='var',
                 colorbar_title='column scaled\nexpression')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/marker_genes.pdf',dpi = 300, bbox_inches='tight')
sc.pl.umap(epi,color=['Pathtype'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Pathtype.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Pathtype.pdf',dpi = 300, bbox_inches='tight')
sc.pl.umap(epi,color=['Dataset'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Dataset.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Dataset.pdf',dpi = 300, bbox_inches='tight')