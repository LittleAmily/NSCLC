from pysam import reference
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
from sympy import use


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
adata = ad1.concatenate(ad2, ad3,ad4, join="inner")
adata.obs['Celltype']
malignant1 = adata[adata.obs['Celltype_fine'] == 'Malignant']
malignant2 = adata[adata.obs['tumor_nontumor_finer'] .isin([ 'Tumor','Tumor (Not Tumor)'])]
malignant3 = adata[adata.obs['cell_type_major'] =='Tumor cells']
malignant = malignant1.concatenate(malignant2,malignant3,join="inner")
malignant.write_h5ad('/home/data/sdzl14/NSCLC/zong/malignant.h5ad')
malignant = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/malignant.h5ad')
malignant = malignant.copy()
# 校正测序深度
sc.pp.normalize_total(malignant, target_sum=1e4)
# 表达值log转化
sc.pp.log1p(malignant)
sc.pp.scale(malignant, max_value=10)
sc.pp.highly_variable_genes(malignant, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes = 2000,batch_key='Dataset')
# 高变基因展示
sc.pl.highly_variable_genes(malignant)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/highly_variable_genes.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/highly_variable_genes.png',dpi = 300, bbox_inches='tight')

sc.pp.pca(malignant,use_highly_variable=True)
sc.pl.pca_variance_ratio(malignant, log=True)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/pca_variance_ratio.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/pca_variance_ratio.png',dpi = 300, bbox_inches='tight')
sc.pp.neighbors(malignant,n_pcs=18)
sc.tl.umap(malignant)
sc.pl.umap(malignant, color=["Dataset"], palette=sc.pl.palettes.vega_20_scanpy, size=50, edgecolor="k", linewidth=0.05, alpha=0.9, show=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/UN_intergrated.pdf")
# List of adata per batch
dataset_cats = malignant.obs['Dataset'].cat.categories.tolist()
malignant_list = [malignant[malignant.obs['Dataset'] == b].copy() for b in dataset_cats]
scanorama.integrate_scanpy(malignant_list)

malignant.obsm["X_scanorama"] = np.zeros((malignant.shape[0], malignant_list[0].obsm["X_scanorama"].shape[1]))
for i,b in enumerate(dataset_cats):
    malignant.obsm["X_scanorama"][malignant.obs['Dataset'] == b] = malignant_list[i].obsm["X_scanorama"]
malignant
sc.pp.neighbors(malignant, use_rep="X_scanorama")
sc.tl.umap(malignant)
sc.set_figure_params(dpi=300,dpi_save=300,figsize=(5,5),vector_friendly=True)
sc.pl.umap(malignant, color=["Dataset"], palette=sc.pl.palettes.vega_20_scanpy, size=50, edgecolor="k", linewidth=0.05, alpha=0.5, show=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/X_scanorama.pdf",dpi = 300,bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/X_scanorama.png",dpi = 300,bbox_inches='tight')
malignant.write_h5ad('/home/data/sdzl14/NSCLC/zong/malignant.h5ad',compression='gzip')
malignant = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/malignant.h5ad')
malignant = malignant.copy()
malignant_scanvi = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/malignant_scANVI.h5ad')
malignant_scanvi = malignant_scanvi.copy()
malignant_scanvi
malignant.obsm['X_scANVI'] = malignant_scanvi.obsm['X_scANVI']
malignant.obsm['X_scVI'] = malignant_scanvi.obsm['X_scVI']
malignant.X.max()
sc.pp.highly_variable_genes(malignant, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes = 2000,batch_key='Dataset')
# 高变基因展示
sc.pl.highly_variable_genes(malignant)
sc.pp.neighbors(malignant,use_rep='X_scANVI')
sc.tl.umap(malignant)
sc.pl.umap(malignant, color=["Dataset"], palette=sc.pl.palettes.vega_20_scanpy, size=50, edgecolor="k", linewidth=0.05, alpha=0.9, show=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Dataset_UMAP.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Dataset_UMAP.png",dpi = 300, bbox_inches='tight')
sc.pl.umap(malignant, color=["Dataset"])
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_UMAP.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_UMAP.png",dpi = 300, bbox_inches='tight')
sc.pl.umap(malignant, color=["Tissue"])
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Tissue_UMAP.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Tissue_UMAP.png",dpi = 300, bbox_inches='tight')
sc.pl.umap(malignant, color=["Origin"])
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Origin_Tissue_UMAP.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Origin_Tissue_UMAP.png",dpi = 300, bbox_inches='tight')
sc.pl.umap(malignant, color=["Drug"])
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Drug_Origin_Tissue_UMAP.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Drug_Origin_Tissue_UMAP.png",dpi = 300, bbox_inches='tight')
sc.pl.umap(malignant, color=["Patient"])
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Patient_Origin_Tissue_UMAP.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Patient_Origin_Tissue_UMAP.png",dpi = 300, bbox_inches='tight')

sc.pl.umap(malignant_scanvi, color=["_scvi_labels"])
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/_scvi_labels.pdf",dpi = 300, bbox_inches='tight')
print(malignant.obs['Origin'].value_counts())
malignant.write_h5ad('/home/data/sdzl14/NSCLC/zong/malignant.h5ad',compression='gzip')
malignant = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/malignant.h5ad')
malignant = malignant.copy()
malignant.obs['Tissue'].cat.categories
adata_spatial = sc.read_h5ad("/home/data/sdzl14/NSCLC/zong/adata_spatial.h5ad")
adata_spatial = adata_spatial.copy()
adata_spatial.var.to_csv("/home/data/sdzl14/NSCLC/zong/spatial/gene_infor.csv")
#读入注释后的基因信息
df = pd.read_csv('/home/data/sdzl14/NSCLC/zong/gene_infor.csv')

#加入基因染色体信息
adata_spatial.var['chromosome'] = df['chr'].to_list()
adata_spatial.var['start'] = df['start'].to_list()
adata_spatial.var['end'] = df['end'].to_list()
#检查基因信息 染色体和基因位置
adata_spatial.var

#去掉NA
adata_spatial = adata[: , adata.var.chromosome.notna()]
#再检查基因信息 染色体和基因位置
adata_spatial.var
import infercnvpy as cnv
##开始跑CNV
cnv.tl.infercnv(
    adata_spatial,
    window_size=250,
    n_jobs= 20
)

##画热图
#设置dpi=600
sc.settings.set_figure_params(dpi=600, facecolor="white")

##画图
#查看不同的分组
#鉴定的细胞类型
cnv.pl.chromosome_heatmap(adata_spatial, groupby="leiden")
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/cnvheatmap_leiden.pdf",dpi = 300, bbox_inches='tight')
#leiden聚类
cnv.tl.pca(adata_spatial)
cnv.pp.neighbors(adata_spatial)
cnv.tl.umap(adata_spatial)
cnv.tl.leiden(adata_spatial, resolution=0.5)
cnv.pl.chromosome_heatmap(adata_spatial, groupby="cnv_leiden")
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/cnvheatmap_cnv_leiden.pdf",dpi = 300, bbox_inches='tight')
cnv.pl.umap(adata_spatial, color="cnv_leiden")
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/cnvheatmap_cnv_leiden_umap.pdf",dpi = 300, bbox_inches='tight')
adata_spatial
cnv.pl.umap(adata_spatial, color="cnv_leiden", legend_loc="on data")
##最后结果在这里
adata.obsm["X_cnv"]
cnv.tl.cnv_score(adata_spatial)
cnv.pl.umap(adata_spatial, color="cnv_score")
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/cnvheatmap_cnv_score.pdf",dpi = 300, bbox_inches='tight')
sc.pl.umap(adata_spatial, color="cnv_score")
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/cnv_score_umap.pdf",dpi = 300, bbox_inches='tight')
adata_spatial.obs['library_id'].cat.categories
map = {'ST_1.1.1':'tumor_middle', 
       'ST_1.2.1':'tumor_edge', 
       'ST_1.2.2':'tumor_edge', 
       'ST_1.3.1':'normal_adjacent', 
       'ST_1.4.1':'normal_distant', 
       'ST_2.1.1':'tumor_middle',
       'ST_2.2.1':'tumor_edge', 
       'ST_2.2.2':'tumor_edge', 
       'ST_3.1.1':'tumor_middle', 
       'ST_3.2.1':'tumor_edge', 
       'ST_3.2.2':'tumor_edge', 
       'ST_3.3.1':'normal_adjacent',
       'ST_5.2.1':'tumor_edge', 
       'ST_5.2.2':'tumor_edge', 
       'ST_5.3.1':'normal_adjacent', 
       'ST_6.2.1':'tumor_edge', 
       'ST_6.3.1':'normal_adjacent', 
       'ST_6.4.1':'normal_distant'
}
adata_spatial.obs['Tissue'] = (
    adata_spatial.obs['library_id']
    .map(map)
    .astype('category')
)
map = {'ST_1.1.1':'P1', 
       'ST_1.2.1':'P1', 
       'ST_1.2.2':'P1', 
       'ST_1.3.1':'P1', 
       'ST_1.4.1':'P1', 
       'ST_2.1.1':'P2',
       'ST_2.2.1':'P2', 
       'ST_2.2.2':'P2', 
       'ST_3.1.1':'P3', 
       'ST_3.2.1':'P3', 
       'ST_3.2.2':'P3', 
       'ST_3.3.1':'P3',
       'ST_5.2.1':'P5', 
       'ST_5.2.2':'P5', 
       'ST_5.3.1':'P5', 
       'ST_6.2.1':'P6', 
       'ST_6.3.1':'P6', 
       'ST_6.4.1':'P6'
}
adata_spatial.obs['Patient'] = (
    adata_spatial.obs['library_id']
    .map(map)
    .astype('category')
)
map = {'ST_1.1.1':'LUSC', 
       'ST_1.2.1':'LUSC', 
       'ST_1.2.2':'LUSC', 
       'ST_1.3.1':'LUSC', 
       'ST_1.4.1':'LUSC', 
       'ST_2.1.1':'LUSC',
       'ST_2.2.1':'LUSC', 
       'ST_2.2.2':'LUSC', 
       'ST_3.1.1':'LUSC', 
       'ST_3.2.1':'LUSC', 
       'ST_3.2.2':'LUSC', 
       'ST_3.3.1':'LUSC',
       'ST_5.2.1':'LUAD', 
       'ST_5.2.2':'LUAD', 
       'ST_5.3.1':'LUAD', 
       'ST_6.2.1':'LUAD', 
       'ST_6.3.1':'LUAD', 
       'ST_6.4.1':'LUAD'
}
adata_spatial.obs['Pathtype'] = (
    adata_spatial.obs['library_id']
    .map(map)
    .astype('category')
)
sc.pl.umap(adata_spatial, color="Tissue")
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/Tissue.pdf",dpi = 300, bbox_inches='tight')
sc.pl.umap(adata_spatial, color="Patient")
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/Patient.pdf",dpi = 300, bbox_inches='tight')
sc.pl.umap(adata_spatial, color="Pathtype")
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/Pathtype.pdf",dpi = 300, bbox_inches='tight')
sc.tl.leiden(adata_spatial,resolution=0.5)
sc.pl.umap(adata_spatial, color="leiden")
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/leiden.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/leiden.png",dpi = 300, bbox_inches='tight')
sc.tl.rank_genes_groups(adata_spatial, "leiden", method="t-test")
sc.pl.rank_genes_groups(adata_spatial, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/rank_genes_groups.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/spatial/rank_genes_groups.png",dpi = 300, bbox_inches='tight')
adata = sc.read_h5ad('/home/data/sdzl14/NSCLC/Peilin_Wang_2025/scRNA_annotated.h5ad')
adata= adata.copy()
adata.obs['Celltype'].cat.categories
adata.var.to_csv("/home/data/sdzl14/NSCLC/zong/spatial/gene_infor.csv")
#读入注释后的基因信息
df = pd.read_csv('/home/data/sdzl14/NSCLC/zong/gene_infor.csv')
#加入基因染色体信息
adata.var['chromosome'] = df['chr'].to_list()
adata.var['start'] = df['start'].to_list()
adata.var['end'] = df['end'].to_list()
#检查基因信息 染色体和基因位置
adata.var

#去掉NA
adata = adata[: , adata.var.chromosome.notna()]
#再检查基因信息 染色体和基因位置
adata.var
cnv.tl.infercnv(
    adata,
    window_size=250,
    n_jobs= 4,
    reference_key='Celltype',
    reference_cat='T'
)
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.umap(adata)
cnv.tl.leiden(adata, resolution=0.5)
cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden")
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/cnvheatmap_cnv_leiden.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/cnvheatmap_cnv_leiden.png",dpi = 300, bbox_inches='tight')
cnv.pl.umap(adata, color="cnv_leiden")
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/cnv_leiden_umap.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/cnv_leiden_umap.png",dpi = 300, bbox_inches='tight')
cnv.tl.cnv_score(adata)
cnv.pl.umap(adata, color="cnv_score")
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/ccnv_score.pdf",dpi = 300, bbox_inches='tight') #这个结果和上面一样
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/ccnv_score.png",dpi = 300, bbox_inches='tight')
sc.pl.violin(adata,keys='cnv_score',
    groupby='cnv_leiden',
    rotation=45,  # 旋转x轴标签角度
    palette='Set2',  # 调色板选择
    figsize=(8, 4),  # 调整图像尺寸
    show=False  # 关闭自动显示（配合保存使用）
)
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/cnv_score_violin.pdf",dpi = 300, bbox_inches='tight') 
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/cnv_score_violin.png",dpi = 300, bbox_inches='tight') 
adata.obs['cnv_score']
cnv.pl.umap(adata, color="Tissue")
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/cnvTissue_.pdf",dpi = 300, bbox_inches='tight')
cnv.pl.umap(adata, color="Patient")
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/cnvPatient_.pdf",dpi = 300, bbox_inches='tight')
cnv.pl.umap(adata, color="Pathtype")
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/cnvPathtype_.pdf",dpi = 300, bbox_inches='tight')
tumor = ['6','7','8','9','11','13']
tumor = adata[adata.obs['cnv_leiden'].isin(tumor)]
tumor1 = malignant[malignant.obs['Dataset'] == 'Peilin_Wang_2025']
tumor2 = malignant[malignant.obs['Patient'] == 'Peilin_Wang_2025_P6']
tumor3 = tumor[tumor.obs['Patient'] == 'Peilin_Wang_2025_P6']
malignant1 = malignant.concatenate(tumor,join="inner")
malignant1.var = {}
malignant1.obs['celltype_fine'].cat.categories
print(tumor.obs['celltype_fine'].value_counts())
sc.tl.rank_genes_groups(adata, 'celltype_fine', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/rank_genes_groups.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/rank_genes_groups.png",dpi = 300, bbox_inches='tight')
adata.obs["cnv_status"] = "normal"
adata.obs.loc[
    adata.obs["cnv_leiden"].isin(['6','7','8','9','11','13']), "cnv_status"
] = "tumor"
cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "tumor", :])
plt.savefig('/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/cnv—leiden分配之后tumor染色体热图.png',dpi=600,bbox_inches='tight')
cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "normal", :])
plt.savefig('/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/cnv—leiden分配之后normal染色体热图.png',dpi=600,bbox_inches='tight')
sc.tl.rank_genes_groups(adata, 'cnv_status', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/rank_genes_groups.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/rank_genes_groups.png",dpi = 300, bbox_inches='tight')
tumor.X.max()
tumor.layers['counts'] = tumor.X
sc.pp.normalize_total(tumor, target_sum=1e4)
sc.pp.log1p(tumor)
sc.pp.highly_variable_genes(tumor, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(tumor)
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/20250417tumor/highly_variable_genes.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/highly_variable_genes.pdf",dpi = 300, bbox_inches='tight')
sc.pp.scale(tumor, max_value=10)
sc.tl.pca(tumor)
sc.pp.neighbors(tumor)
sc.tl.umap(tumor)
sc.tl.leiden(tumor, resolution=0.1)
sc.pl.umap(tumor, color="leiden")
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/20250417tumor/umap.png",dpi = 300, bbox_inches='tight')
sc.tl.rank_genes_groups(tumor, 'leiden', method='t-test')
sc.pl.rank_genes_groups(tumor, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/20250417tumor/rank_genes_groups.png",dpi = 300, bbox_inches='tight')
print(tumor.obs['celltype_fine'].value_counts())
print(tumor3.obs['celltype_fine'].value_counts())
sc.pl.umap(tumor, color="celltype_main")
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/20250417tumor/umap_celltype_maine.png",dpi = 300, bbox_inches='tight')
sc.pl.umap(tumor, color="Patient")
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/20250417tumor/umap_Patient.png",dpi = 300, bbox_inches='tight')
sc.pl.umap(tumor, color="Tissue")
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/20250417tumor/umap_Tissue.png",dpi = 300, bbox_inches='tight')
sc.pl.umap(tumor3, color="celltype_main")
plt.savefig("/home/data/sdzl14/NSCLC/Peilin_Wang_2025/pdf/20250417tumor/umap3_celltype_main.png",dpi = 300, bbox_inches='tight')
malignant.obs['Tissue'].cat.categories
map = {'effusion':'effusion', 
       'normal_adjacent':'normal_adjacent', 
       'normal_adjancent':'normal_adjacent', 
       'normal_distant':'normal_distant',
       'tumor_edge':'tumor_edge', 
       'tumor_metastasis':'tumor_metastasis', 
       'tumor_metastatis':'tumor_metastasis', 
       'tumor_middle':'tumor_middle',
       'tumor_primary':'tumor_primary'}
malignant.obs['Tissue'] = (
    malignant.obs['Tissue']
    .map(map)
    .astype('category')
)
sc.pl.umap(malignant, color="Tissue")
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Tissue_UMAP.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Tissue_UMAP.png",dpi = 300, bbox_inches='tight')
sc.tl.leiden(malignant, resolution=0.1)
sc.pl.umap(malignant, color="leiden")
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_leiden_UMAP.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_leiden_UMAP.png",dpi = 300, bbox_inches='tight')
LUSC = malignant[malignant.obs['Pathtype'] == 'LUSC']
sc.tl.rank_genes_groups(LUSC, 'Tissue', method='t-test')
sc.pl.rank_genes_groups(LUSC, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_intergrated_rank_genes_groups.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_intergrated_rank_genes_groups.png",dpi = 300, bbox_inches='tight')
LUSC.obs['Tissue'].value_counts()
sc.pp.highly_variable_genes(LUSC, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(LUSC)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_highly_variable_genes.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_highly_variable_genes.pdf",dpi = 300, bbox_inches='tight')
sc.pp.scale(LUSC, max_value=10)
sc.tl.pca(LUSC)
sc.pp.neighbors(LUSC,use_rep='X_scANVI')
sc.tl.umap(LUSC)
sc.tl.leiden(LUSC, resolution=0.5)
sc.pl.umap(LUSC, color="leiden")
map = {
    '0':'0',
    '1':'1',
    '2':'2',
    '3':'3',
    '4':'4',
    '5':'5',
    '6':'6',
    '7':'7',
    '8':'4',
    '9':'4',
    '10':'0',
    '11':'1',
    '12':'12',
}
LUSC.obs['leiden'] = (
    LUSC.obs['leiden']
    .map(map)
    .astype('category')
)
LUSC.obs['Tissue'] = (
    LUSC.obs['Tissue']
    .map(map)
    .astype('category')
)
lusc = ['0', '1', '2', '3', '4', '5', '6', '7']
LUSC = LUSC[LUSC.obs['leiden'].isin(lusc)]
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap.pdf",dpi = 300, bbox_inches='tight')
sc.pl.umap(LUSC, color="Tissue")
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Tissue.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Tissue.pdf",dpi = 300, bbox_inches='tight')
sc.pl.umap(LUSC, color="Patient")
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Patient.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Patient.pdf",dpi = 300, bbox_inches='tight')
sc.pl.umap(LUSC, color="Origin")
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Origin.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Origin.pdf",dpi = 300, bbox_inches='tight')
sc.pl.umap(LUSC, color="Timepoint")
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Timepoint.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Timepoint.pdf",dpi = 300, bbox_inches='tight')
import pandas as pd
from anndata import AnnData

# 假设你的adata是AnnData对象
def calculate_cell_percentage(adata: AnnData, 
                             leiden_col: str = "leiden",  # 根据实际情况修改列名
                             tissue_col: str = "Tissue",
                             output_path: str = "/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_cell_counts_Originpercentage.csv"):
    # 生成统计计数
    cell_counts = (
        adata.obs.groupby([tissue_col, leiden_col])
        .size()
        .reset_index(name='n')
    )
    
    # 计算各组织总细胞数
    total_cells = (
        cell_counts.groupby(tissue_col)
        ['n'].sum()
        .reset_index(name='total')
    )
    
    # 合并统计量并计算百分比
    merged_df = pd.merge(
        cell_counts, 
        total_cells, 
        on=tissue_col,
        how='left'
    )
    
    merged_df['percentage'] = merged_df['n'] / merged_df['total'] * 100
    
    # 保存结果
    merged_df.to_csv(output_path, index=False)
    return merged_df

df_result = calculate_cell_percentage(LUSC, leiden_col='leiden', output_path='/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_Origincell_counts_percentage.csv')
# 验证每个组织的百分比总和应为100%
check = df_result.groupby('Tissue')['percentage'].sum()
print("\n组织类型百分比验证:")
print(check.round(2))
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_leiden_distribution(df: pd.DataFrame, 
                            output_path: str = "leiden_distribution.pdf",
                            palette: list = None,
                            figsize: tuple = (10, 6)):
    """
    参数说明：
    df: 必须包含的列 ['Tissue', 'leiden', 'n', 'percentage']
    palette: 颜色列表，长度需与Tissue类别数一致
    """
    
    # ================ 新增排序预处理 ================
    # 定义强制排序逻辑
    tissue_order = ['tumor_middle', 'tumor_edge', 
                   'normal_adjacent', 'tumor_metastasis']
    
    # 转换为分类类型确保排序
    df['Tissue'] = pd.Categorical(
        df['Tissue'], 
        categories=tissue_order,
        ordered=True
    )
    
    # 过滤不在排序列表中的异常数据（可选）
    df = df[df['Tissue'].isin(tissue_order)].copy()
    
    # ================ 绘图设置 ================
    sns.set_theme(style="white", font_scale=1.2)
    plt.rcParams.update({
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10
    })

    # ================ 左侧百分比堆积图 ================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, 
                                  gridspec_kw={'width_ratios': [2, 1]})

    # 数据透视时保持排序
    pivot_df = df.pivot_table(
        index='Tissue',
        columns='leiden',
        values='percentage',
        aggfunc='sum'
    ).fillna(0).reindex(tissue_order)  # 关键排序操作

    # 绘制堆叠条形图
    pivot_df.plot.barh(
        stacked=True,
        ax=ax1,
        color=palette,
        edgecolor='black',
        linewidth=0.5,
        width=0.8
    )
    
    # 调整坐标轴方向
    ax1.invert_yaxis()  # 使tumor_middle显示在最上方
    ax1.set(
        ylabel='Tissue Category',  # 修改标签更准确
        xlabel='Percentage (%)',
        xlim=(0, 100)
    )
    # ================ 右侧细胞数量图 ================
    # 按排序后的tissue顺序计算总数
    total_df = df.groupby('leiden', observed=True)['n'].sum().reset_index()
    
    sns.barplot(
        data=total_df,
        y='leiden',
        x='n',
        ax=ax2,
        color='#008B8B',
        edgecolor='black',
        linewidth=0.5,
        order=sorted(df['leiden'].unique())  # 按leiden数字排序
    )
    # 样式调整
    ax2.set(
        xlabel='Cell Count', 
        ylabel='',
        yticklabels=[]
    )
    ax2.grid(False)
    
    # ====================== 整体布局调整 ======================
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()

# 使用示例（需要替换实际参数）：
# 假设df_result是已计算好的数据框
tissue_palette = ['#fcf1f0', '#fccccb','#bdb5e1','#b0d992','#f9d580','#85a2d2','#e3716e','#eca680','#7ac7e2','#f7df87','#54beaa','#2983b1','#882E71']  # 颜色数量需与Tissue类别数一致
plot_leiden_distribution(df_result, 
                        palette=tissue_palette,
                        output_path='/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_leiden_distribution--.pdf',
                        figsize=(9, 6))
import pandas as pd
from anndata import AnnData



# 验证每个组织的百分比总和应为100%
check = df_result.groupby('leiden')['percentage'].sum()
print("\nleiden类型百分比验证:")
print(check.round(2))
plot_leiden_distribution(
    df_result,

    output_path='/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_tissue_leiden_distribution.pdf',
    figsize=(10, 6)
)

sc.tl.rank_genes_groups(LUSC, "leiden", method="t-test")
sc.pl.rank_genes_groups(LUSC, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_rank_genes_groups.png", dpi=300, bbox_inches='tight')
LUSC.obs['Tissue'].cat.categories
LUSC_ = LUSC[LUSC.obs['Tissue'].isin(['normal_adjacent', 'tumor_edge', 'tumor_metastasis',
       'tumor_middle'])]
sc.tl.rank_genes_groups(LUSC_, "leiden", method="t-test")
sc.pl.rank_genes_groups(LUSC_, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_rank_genes_groups.png", dpi=300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_rank_genes_groups.pdf", dpi=300, bbox_inches='tight')
sc.tl.rank_genes_groups(LUSC_, "Tissue", method="t-test")
sc.pl.rank_genes_groups(LUSC_, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_intergrated_rank_genes_groups.pdf", dpi=300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_intergrated_rank_genes_groups.png", dpi=300, bbox_inches='tight')
def plot_leiden2_distribution(df_result, output_path, figsize=(6,6)):
    # 过滤数据：仅保留leiden=2的数据
    df_leiden2 = df_result[df_result['leiden'] == '2'].copy()
    
    # 设置自定义排序逻辑
    tissue_order = ['tumor_middle', 'tumor_edge', 
                   'normal_adjacent', 'tumor_metastasis']
    
    # 转换为分类类型确保排序
    df_leiden2['Tissue'] = pd.Categorical(
        df_leiden2['Tissue'], 
        categories=tissue_order,
        ordered=True
    )
    
    # 创建画布
    plt.figure(figsize=figsize, dpi=300)
    
    # 绘制折线图
    sns.lineplot(
        data=df_leiden2.sort_values('Tissue'),  # 强制按自定义顺序排序
        x='Tissue',
        y='percentage',
        marker='o',
        markersize=8,
        linewidth=2.5,
        color='#2983b1'  # 使用蓝色系颜色
    )
    
    # 设置坐标轴标签
    plt.xlabel('Tissue Category', fontsize=12)
    plt.ylabel('Percentage of Leiden 2 (%)', fontsize=12)
    plt.title('Leiden Cluster 2 Proportion Across Tissues', fontsize=14, pad=20)
    
    # 调整x轴标签
    plt.xticks(
        rotation=45,
        ha='right',
        fontsize=11,
        rotation_mode='anchor'
    )
    
    # 设置y轴范围
    plt.ylim(0, df_leiden2['percentage'].max()*1.1)
    
    # 添加数值标签
    for x, y in zip(df_leiden2['Tissue'], df_leiden2['percentage']):
        plt.text(
            x=x, 
            y=y+0.5,  # 向上偏移0.5%
            s=f"{y:.1f}%",
            ha='center',
            va='bottom',
            fontsize=10
        )
    
    # 优化布局
    plt.tight_layout()
    
    # 保存输出
    plt.savefig(output_path, bbox_inches='tight')
    plt.savefig(output_path.replace('.pdf','.png'), bbox_inches='tight')
    plt.close()

# 调用函数
plot_leiden2_distribution(
    df_result,
    output_path='/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_leiden2_distribution.pdf',
    figsize=(5, 5)
)
import palantir
## (1)Run diffusion maps
pca_projections = pd.DataFrame(LUSC_.obsm['X_pca'], index=LUSC_.obs_names)
umap = pd.DataFrame(LUSC_.obsm['X_umap'], index=LUSC_.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=5)

ms_data = palantir.utils.determine_multiscale_space(dm_res)
LUSC_.layers['MAGIC_imputed_data'] = palantir.utils.run_magic_imputation(LUSC_, dm_res)

#基因表达量可视化
sc.pl.embedding(LUSC_, basis='umap', layer='MAGIC_imputed_data',
 color=['S100A9', 'S100A8', 'S100A4', 'S100A2','S100A10'])
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_magic_基因表达量可视化.png", dpi=300, bbox_inches='tight')

# 可视化
sc.pl.embedding(
adata,
color=["ct_pseudotime", "Stage"],
basis="force_directed",
color_map="gnuplot2",
)
# 设置目标基因和组织类型
target_genes = ['SEC61G','NME2','IGLC2','APOC1']
tissue_order = ['tumor_middle', 'tumor_edge', 'normal_adjacent', 'tumor_metastasis']

# 提取表达矩阵和元数据
expression_df = pd.DataFrame(
    LUSC_[:, target_genes].X.toarray(),  # 假设X是稀疏矩阵
    columns=target_genes,
    index=LUSC_.obs_names
)
meta_df = LUSC_.obs[['Tissue']].copy()

# 合并数据
combined_df = pd.concat([expression_df, meta_df], axis=1)

# 计算平均表达量（按组织分组）
mean_expression = (
    combined_df.groupby('Tissue')[target_genes]
    .mean()
    .reset_index()
    .melt(id_vars='Tissue', var_name='Gene', value_name='Expression')
)

# 强制排序
mean_expression['Tissue'] = pd.Categorical(
    mean_expression['Tissue'], 
    categories=tissue_order,
    ordered=True
)

# 可视化设置
plt.figure(figsize=(6, 6), dpi=300)
sns.set_style("whitegrid", {'grid.linestyle': '--'})

# 设置颜色方案
gene_palette = [ '#fccccb','#bdb5e1','#b0d992','#f9d580']

# 修改绘图部分代码
ax = sns.lineplot(
    data=mean_expression.sort_values('Tissue'),
    x='Tissue',
    y='Expression',
    hue='Gene',
    palette=gene_palette,  # 使用自定义颜色
    marker='o',
    markersize=10,  # 增大标记尺寸
    linewidth=3,     # 加粗线条
    err_style='bars',
    errorbar=('se', 1.96)
)

# 优化图例显示
plt.legend(
    title='Gene',
    frameon=True,
    facecolor='white',
    bbox_to_anchor=(1.25, 0.5),
    title_fontsize=12,
    fontsize=11
)

# 样式优化
plt.title("S100A Family Gene Expression Across Tissue Types", 
         fontsize=14, pad=20)
plt.xlabel("Tissue Category", fontsize=12)
plt.ylabel("Log Normalized Expression", fontsize=12)
plt.xticks(rotation=45, ha='right', fontsize=11)
plt.legend(title='Gene', frameon=False, bbox_to_anchor=(1, 0.5))

# 保存输出
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_expression.pdf', 
           bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_expression.png', 
           dpi=300, bbox_inches='tight')
plt.close()
LUAD = malignant[malignant.obs['Pathtype'] == 'LUAD']
malignant = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/malignant.h5ad')
malignant = malignant.copy()
malignant.obs['Tissue'].cat.categories
mali = ['effusion', 'normal_adjacent', 'normal_adjancent', 'normal_distant',
       'tumor_edge', 'tumor_metastasis', 'tumor_metastatis', 'tumor_middle']

malignant = malignant[malignant.obs['Tissue'].isin(mali)]
map = {
    'effusion':'effusion', 
    'normal_adjacent':'normal_adjacent', 
    'normal_adjancent':'normal_adjacent',
    'normal_distant':'normal_distant',
    'tumor_edge':'tumor_edge', 
    'tumor_metastasis':'tumor_metastasis', 
    'tumor_metastatis':'tumor_metastasis', 
    'tumor_middle':'tumor_middle'
}

malignant.obs['Tissue'] = (
    malignant.obs['Tissue']
    .map(map)
    .astype('category')
)
malignant.X.max()

sc.pp.highly_variable_genes(
    malignant,
    flavor="seurat",
    n_top_genes=2000,
    inplace=True
        )
sc.pp.pca(malignant,use_highly_variable=True)
sc.pl.pca_variance_ratio(malignant, log=True)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/pca_variance_ratio.pdf',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/pca_variance_ratio.png',dpi = 300, bbox_inches='tight')
sc.pp.neighbors(malignant,use_rep='X_scANVI')
sc.tl.umap(malignant)
sc.pl.umap(malignant, color=["Dataset"],show=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_UMAP.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_UMAP.png",dpi = 300, bbox_inches='tight')
malignant.obs['Dataset'].cat.categories
malignant = malignant[malignant.obs['Dataset'].isin(['Goveia_Carmeliet_2020', 'He_Fan_2021', 'Kim_Lee_2020',
       'Lambrechts_Thienpont_2018_6149v1', 'Lambrechts_Thienpont_2018_6149v2',
       'Lambrechts_Thienpont_2018_6653', 'Laughney_Massague_2020',
       'Maynard_Bivona_2020', 'Peilin_Wang_2025', 
       'Travaglini_Krasnow_2020', 'UKIM-V', 'Vieira_Teichmann_2019'])]
sc.pp.neighbors(malignant,use_rep='X_scANVI')
sc.tl.umap(malignant)
sc.pl.umap(malignant, color=["Dataset"],show=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_UMAP.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_UMAP.png",dpi = 300, bbox_inches='tight')
sc.tl.leiden(malignant, resolution=0.2)
sc.pl.umap(malignant, color=["leiden"],show=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_leiden.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_leiden.png",dpi = 300, bbox_inches='tight')
# 假设你的adata是AnnData对象
def calculate_cell_percentage(adata: AnnData, 
                             leiden_col: str = "leiden",  # 根据实际情况修改列名
                             tissue_col: str = "Tissue",
                             output_path: str = "/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_cell_counts_Originpercentage.csv"):
    # 生成统计计数
    cell_counts = (
        adata.obs.groupby([tissue_col, leiden_col])
        .size()
        .reset_index(name='n')
    )
    
    # 计算各组织总细胞数
    total_cells = (
        cell_counts.groupby(tissue_col)
        ['n'].sum()
        .reset_index(name='total')
    )
    
    # 合并统计量并计算百分比
    merged_df = pd.merge(
        cell_counts, 
        total_cells, 
        on=tissue_col,
        how='left'
    )
    
    merged_df['percentage'] = merged_df['n'] / merged_df['total'] * 100
    
    # 保存结果
    merged_df.to_csv(output_path, index=False)
    return merged_df

df_result = calculate_cell_percentage(malignant, leiden_col='leiden', output_path='/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/inte_Origincell_counts_percentage.csv')
# 验证每个组织的百分比总和应为100%
check = df_result.groupby('Tissue')['percentage'].sum()
print("\n组织类型百分比验证:")
print(check.round(2))
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_leiden_distribution(df: pd.DataFrame, 
                            output_path: str = "leiden_distribution.pdf",
                            palette: list = None,
                            figsize: tuple = (10, 6)):
    """
    参数说明：
    df: 必须包含的列 ['Tissue', 'leiden', 'n', 'percentage']
    palette: 颜色列表，长度需与Tissue类别数一致
    """
    
    # ================ 新增排序预处理 ================
    # 定义强制排序逻辑
    tissue_order = ['tumor_middle', 'tumor_edge', 
                   'normal_adjacent', 'tumor_metastasis']
    
    # 转换为分类类型确保排序
    df['Tissue'] = pd.Categorical(
        df['Tissue'], 
        categories=tissue_order,
        ordered=True
    )
    
    # 过滤不在排序列表中的异常数据（可选）
    df = df[df['Tissue'].isin(tissue_order)].copy()
    
    # ================ 绘图设置 ================
    sns.set_theme(style="white", font_scale=1.2)
    plt.rcParams.update({
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10
    })

    # ================ 左侧百分比堆积图 ================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, 
                                  gridspec_kw={'width_ratios': [2, 1]})

    # 数据透视时保持排序
    pivot_df = df.pivot_table(
        index='Tissue',
        columns='leiden',
        values='percentage',
        aggfunc='sum'
    ).fillna(0).reindex(tissue_order)  # 关键排序操作

    # 绘制堆叠条形图
    pivot_df.plot.barh(
        stacked=True,
        ax=ax1,
        color=palette,
        edgecolor='black',
        linewidth=0.5,
        width=0.8
    )
    
    # 调整坐标轴方向
    ax1.invert_yaxis()  # 使tumor_middle显示在最上方
    ax1.set(
        ylabel='Tissue Category',  # 修改标签更准确
        xlabel='Percentage (%)',
        xlim=(0, 100)
    )
    # ================ 右侧细胞数量图 ================
    # 按排序后的tissue顺序计算总数
    total_df = df.groupby('leiden', observed=True)['n'].sum().reset_index()
    
    sns.barplot(
        data=total_df,
        y='leiden',
        x='n',
        ax=ax2,
        color='#008B8B',
        edgecolor='black',
        linewidth=0.5,
        order=sorted(df['leiden'].unique())  # 按leiden数字排序
    )
    # 样式调整
    ax2.set(
        xlabel='Cell Count', 
        ylabel='',
        yticklabels=[]
    )
    ax2.grid(False)
    
    # ====================== 整体布局调整 ======================
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()

# 使用示例（需要替换实际参数）：
# 假设df_result是已计算好的数据框
tissue_palette = ['#fcf1f0', '#fccccb','#bdb5e1','#b0d992','#f9d580','#85a2d2','#e3716e','#eca680','#7ac7e2','#f7df87','#54beaa','#2983b1','#882E71']  # 颜色数量需与Tissue类别数一致
plot_leiden_distribution(df_result, 
                        palette=tissue_palette,
                        output_path='/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_leiden_distribution--.pdf',
                        figsize=(9, 6))
sc.tl.rank_genes_groups(malignant, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(malignant, n_genes=25, sharey=False, show=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_rank_genes_groups.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_rank_genes_groups.pdf",dpi = 300, bbox_inches='tight')
sc.tl.rank_genes_groups(malignant, 'leiden', method='t-test')
sc.pl.rank_genes_groups(malignant, n_genes=25, sharey=False, show=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_rank_genes_groups-.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_rank_genes_group-.pdf",dpi = 300, bbox_inches='tight')
sc.tl.rank_genes_groups(malignant, 'Tissue', method='wilcoxon')
sc.pl.rank_genes_groups(malignant, n_genes=25, sharey=False, show=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Tissuerank_genes_groups.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Tissuerank_genes_groups.pdf",dpi = 300, bbox_inches='tight')
sc.tl.rank_genes_groups(malignant, 'Tissue', method='t-test')
sc.pl.rank_genes_groups(malignant, n_genes=25, sharey=False, show=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Tissuerank_genes_groups-.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/intergrated_Tissuerank_genes_groups-.pdf",dpi = 300, bbox_inches='tight')
# 设置目标基因和组织类型
target_genes = ['NME2','IGLC2','ALDOA','CLU']

tissue_order = ['tumor_middle', 'tumor_edge', 'normal_adjacent', 'tumor_metastasis']

# 提取表达矩阵和元数据
expression_df = pd.DataFrame(
    LUAD[:, target_genes].X.toarray(),  # 假设X是稀疏矩阵
    columns=target_genes,
    index=LUAD.obs_names
)
meta_df = LUAD.obs[['Tissue']].copy()

# 合并数据
combined_df = pd.concat([expression_df, meta_df], axis=1)

# 计算平均表达量（按组织分组）
mean_expression = (
    combined_df.groupby('Tissue')[target_genes]
    .mean()
    .reset_index()
    .melt(id_vars='Tissue', var_name='Gene', value_name='Expression')
)

# 强制排序
mean_expression['Tissue'] = pd.Categorical(
    mean_expression['Tissue'], 
    categories=tissue_order,
    ordered=True
)

# 可视化设置
plt.figure(figsize=(6, 6), dpi=300)
sns.set_style("whitegrid", {'grid.linestyle': '--'})

# 设置颜色方案
gene_palette = ['#fcf1f0', '#fccccb','#bdb5e1','#b0d992','#f9d580','#85a2d2','#e3716e','#eca680','#7ac7e2','#f7df87','#54beaa','#2983b1','#882E71'] 
# 修改绘图部分代码
ax = sns.lineplot(
    data=mean_expression.sort_values('Tissue'),
    x='Tissue',
    y='Expression',
    hue='Gene',
    palette=gene_palette,  # 使用自定义颜色
    marker='o',
    markersize=10,  # 增大标记尺寸
    linewidth=3,     # 加粗线条
    err_style='bars',
    errorbar=('se', 1.96)
)

# 优化图例显示
plt.legend(
    title='Gene',
    frameon=True,
    facecolor='white',
    bbox_to_anchor=(1.25, 0.5),
    title_fontsize=12,
    fontsize=11
)

# 样式优化
plt.title(" Gene Expression Across Tissue Types", 
         fontsize=14, pad=20)
plt.xlabel("Tissue Category", fontsize=12)
plt.ylabel("Log Normalized Expression", fontsize=12)
plt.xticks(rotation=45, ha='right', fontsize=11)
plt.legend(title='Gene', frameon=False, bbox_to_anchor=(1, 0.5))

# 保存输出
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/integratedLUAD_expression.pdf', 
           bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/integratedLUAD_expression.png', 
           dpi=300, bbox_inches='tight')
plt.close()
LUSC = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/LUSC.h5ad')
LUSC = LUSC.copy()
LUSC
LUAD = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/LUAD_qc.h5ad')
LUAD = LUAD.copy()
LUAD
adata_spatial = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/LUAD_spatial.h5ad')
malignant.write_h5ad('/home/data/sdzl14/NSCLC/zong/malignant_integrated.h5ad')   