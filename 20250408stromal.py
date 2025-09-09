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
ad2
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
ad1.obs['celltype_fine'].cat.categories
stromal1 = ['Fibroblasts','Endothelial cells']
stromal1 = ad1[ad1.obs['celltype_fine'].isin(stromal1)]
ad2.obs['minor'].cat.categories
stromal2 = ['ADH1B+ CAF', 'AM', 'Artery', 'BCHE+ SMC', 'COL11A1+ CAF',
       'CPE+ Venule', 'Capillary', 'Lymphatic EC', 'MYH11+ Pericyte','Pericyte','SELE+ Venule', 'SMC', 'Tip','Venule' ]
stromal2 = ad2[ad2.obs['minor'].isin(stromal2)]
ad3.obs['cell_type'].cat.categories
stromal3 = ['Endothelial cell arterial', 'Endothelial cell capillary','Mesothelial',
       'Endothelial cell lymphatic', 'Endothelial cell venous','Fibroblast adventitial', 'Fibroblast alveolar','Pericyte',
       'Fibroblast peribronchial','Smooth muscle cell','stromal dividing']
stromal3 = ad3[ad3.obs['cell_type'].isin(stromal3)]
ad4.obs['tumor_nontumor_finer'].cat.categories
stromal4 = ['Endothelial','Other','Fibroblasts']
stromal4 = ad4[ad4.obs['tumor_nontumor_finer'].isin(stromal4)]
stromal = stromal1.concatenate(stromal2,stromal3,stromal4,join="inner")
keep_obs_columns = ['Sample', 'Patient', 'Celltype', 'Dataset', 'Platform', 
                    'Pathtype', 'Drug', 'Timepoint', 'Tissue', 'Origin','tumor_nontumor_finer','minor','cell_type','celltype_fine']
stromal.obs = stromal.obs[keep_obs_columns].copy()

# 清理 var 的元数据（保留基因名，删除所有列）
stromal.var = pd.DataFrame(index=stromal.var.index)  # 仅保留基因名，删除所有列

# 删除 obsm 中的嵌入数据
stromal.obsm = {}

# 删除 layers 中的额外数据层
stromal.layers = {}
stromal.uns = {}
stromal.obsp = {}
stromal.write_h5ad("/home/data/sdzl14/NSCLC/zong/stromal.h5ad")
