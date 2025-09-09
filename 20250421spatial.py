import os
os.environ["OMP_NUM_THREADS"] = "4"  # 添加在文件最开
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import AnnData
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns


import scanpy  as sc
import os
results_folder = '/home/data/sdzl14/NSCLC/zong/spatial/'
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f"{results_folder}/reference_signatures"
run_name = f"{results_folder}/cell2location_map"
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
adata_vis.obs[adata_vis.uns["mod"]["factor_names"]] = adata_vis.obsm["q05_cell_abundance_w_sf"]
adata_spatial = adata_vis.copy()
adata_spatial.layers['counts'] = adata_spatial.X
sc.pp.normalize_total(adata_spatial, inplace=True)
sc.pp.log1p(adata_spatial)
sc.pp.highly_variable_genes(adata_spatial, flavor="seurat", n_top_genes=2000, inplace=True)
sc.tl.pca(adata_spatial)

import scanpy.external as sce
sce.pp.harmony_integrate(adata_spatial, 'library_id')
sc.pp.neighbors(adata_spatial,use_rep='X_pca_harmony')
sc.tl.umap(adata_spatial)
sc.tl.leiden(adata_spatial,resolution=0.8)
sc.pl.umap(adata_spatial,color=['leiden'],edges=False, legend_loc='on data')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/leiden.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/leiden.png",dpi = 300, bbox_inches='tight')
# 定义要绘制的颜色列表
from turtle import width

from shapely import box


colors = ['AT1', 'AT2', 'CD4+T', 'CD8+T', 'Cancer_P1/2/3', 'Cancer_P4', 'Cancer_P5', 
          'Ciliated cells', 'Club cells', 'Cycling cells', 'DC', 'Endothelial cells', 
          'Epithelial cells', 'Fibroblasts', 'Macrophage', 'Mast cells', 'Memory B cells', 
          'Monocytes', 'NK cells', 'Naive B cells', 'Neutrophils', 'Plasma cells', 'pDC']

# 计算需要的子图数量和行数
num_plots = len(colors)
ncols = 4
nrows = (num_plots + ncols - 1) // ncols  # 向上取整

# 创建子图
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 5 * nrows))

# 确保axs是一个二维数组
if nrows == 1:
    axs = axs.reshape(1, -1)

# 绘制每个UMAP图
for i, color in enumerate(colors):
    row = i // ncols
    col = i % ncols
    sc.pl.umap(adata_spatial, color=color, edges=False, ax=axs[row, col], show=False)

# 删除多余的子图
for i in range(num_plots, nrows * ncols):
    row = i // ncols
    col = i % ncols
    fig.delaxes(axs[row, col])

# 调整子图之间的间距
plt.tight_layout()

# 显示图表
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_cell_score.png",dpi = 300,bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_cell_score.pdf",dpi = 300,bbox_inches='tight')
# 定义细胞类型列表（与obs中的列名严格对应）
cell_type_columns = [
    'AT1', 'AT2', 'CD4+T', 'CD8+T', 'Cancer_P1/2/3', 'Cancer_P4', 'Cancer_P5',
    'Ciliated cells', 'Club cells', 'Cycling cells', 'DC', 'Endothelial cells',
    'Epithelial cells', 'Fibroblasts', 'Macrophage', 'Mast cells', 'Memory B cells',
    'Monocytes', 'NK cells', 'Naive B cells', 'Neutrophils', 'Plasma cells', 'pDC'
]

# 核心操作：对每个spot选择评分最高的细胞类型
adata_spatial.obs['spatial_celltype_fine'] = adata_spatial.obs[cell_type_columns].idxmax(axis=1)
sc.pl.umap(adata_spatial,color=['spatial_celltype_fine'],edges=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/spatial_celltype_fine.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/spatial_celltype_fine.png",dpi = 300, bbox_inches='tight')
adata_spatial.obs['spatial_celltype'].cat.categories
map = {
    'AT1': 'Normal_Epi',
    'AT2': 'Normal_Epi',
    'CD4+T': 'NK&T',
    'CD8+T': 'NK&T',
    'Cancer_P1/2/3': 'LUSC_Malignant',
    'Cancer_P4': 'LUSC_Malignant',
    'Cancer_P5': 'LUAD_Malignant',
    'Ciliated cells': 'Normal_Epi',
    'Club cells': 'Normal_Epi',
    'Cycling cells': 'Normal_Epi',
    'Endothelial cells':'Stromal_Cell',
    'Fibroblasts': 'Stromal_Cell',
    'Macrophage': 'Macrophage',
    'Memory B cells':'B',
    'NK cells':'NK&T',
    'Plasma cells':'B'
}
adata_spatial.obs['spatial_celltype_fine'] = adata_spatial.obs['spatial_celltype'].map(map)
sc.pl.umap(adata_spatial,color=['spatial_celltype'],edges=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/spatial_celltype.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/spatial_celltype.png",dpi = 300, bbox_inches='tight')
map = {
    '0':'Normal_Epi',
    '1':'LUSC_Malignant',
    '2':'LUAD_Malignant',
    '3':'LUSC_Malignant',
    '4':'LUSC_Malignant',
    '5':'LUSC_Malignant',
    '6':'LUAD_Malignant',
    '7':'B',
    '8':'Stromal_Cell',
    '9':'Macrophage',
    '10':'Ciliated Cell',
    '11':'Macrophage',
    '12':'B',
    '13':'LUSC_Malignant'
    
}
adata_spatial.obs['spatial_celltype_coarse'] = adata_spatial.obs['leiden'].map(map)
sc.pl.umap(adata_spatial,color=['spatial_celltype_coarse'],edges=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/spatial_celltype_coarse.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/spatial_celltype_coarse.png",dpi = 300, bbox_inches='tight')
adata_spatial.obs['library_id'].cat.categories
library_id = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1', 'ST_2.1.1',
       'ST_2.2.1', 'ST_2.2.2', 'ST_3.1.1', 'ST_3.2.1', 'ST_3.2.2', 'ST_3.3.1',
       'ST_5.2.1', 'ST_5.2.2', 'ST_5.3.1', 'ST_6.2.1', 'ST_6.3.1', 'ST_6.4.1']
map = {
    'ST_1.1.1':'P1',
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
adata_spatial.obs['Patient'] = adata_spatial.obs['library_id'].map(map)
sc.pl.umap(adata_spatial,color=['Patient'],edges=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/Patient.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/Patient.png",dpi = 300, bbox_inches='tight')
map = {
    'ST_1.1.1':'tumor_middle',
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
adata_spatial.obs['Tissue'] = adata_spatial.obs['library_id'].map(map)
sc.pl.umap(adata_spatial,color=['Tissue'],edges=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/Tissue.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/Tissue.png",dpi = 300, bbox_inches='tight')
sc.pl.umap(adata_spatial,color=['NME2','IGLC2','CLU','ALDOA'],edges=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/GeneExpression.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/GeneExpression.png",dpi = 300, bbox_inches='tight')
# 定位需要修改的条件
# 添加新类别 'Cancer_P6'
adata_spatial.obs['spatial_celltype_fine'] = adata_spatial.obs['spatial_celltype_fine'].cat.add_categories(['Cancer_P6'])
condition = (adata_spatial.obs['Patient'] == 'P6') & \
            (adata_spatial.obs['spatial_celltype_fine'] == 'Cancer_P1/2/3')

# 执行修改
adata_spatial.obs.loc[condition, 'spatial_celltype_fine'] = 'Cancer_P6'
library_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1', 'ST_2.1.1',
       'ST_2.2.1', 'ST_2.2.2', 'ST_3.1.1', 'ST_3.2.1', 'ST_3.2.2', 'ST_3.3.1',
       'ST_5.2.1', 'ST_5.2.2', 'ST_5.3.1', 'ST_6.2.1', 'ST_6.3.1', 'ST_6.4.1']

for lib_id in library_ids:
    # 过滤当前样本的数据
    adata_subset = adata_spatial[adata_spatial.obs['library_id'] == lib_id, :].copy()
    
    # 关键修改：将 library_id 参数替换为动态变量 lib_id
    sc.pl.spatial(adata_subset, 
                 img_key="hires",
                 color=['YBX1'],
                 edges=False,
                 library_id=lib_id,  # 动态传递当前样本ID
                 title=f"YBX1 in {lib_id}",
                 show=False)
    
    # 保存文件
    plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/Gene_Expression_YBX1_{lib_id}.pdf", 
            dpi=300, bbox_inches='tight')
    plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/Gene_Expression_YBX1_{lib_id}.png",
            dpi=300, bbox_inches='tight')
    plt.close()
sc.pl.spatial(adata_spatial,img_key="hires",color=['CLU'],edges=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/GeneExpression_CLU.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/GeneExpression_CLU.png",dpi = 300, bbox_inches='tight')
sc.pl.spatial(adata_spatial,img_key="hires",color=['ALDOA'],edges=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/GeneExpression_ALDOA.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/GeneExpression_ALDOA.png",dpi = 300, bbox_inches='tight')
sc.pl.spatial(adata_spatial,img_key="hires",color=['NME2'],edges=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/GeneExpression_NME2.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/GeneExpression_NME2.png",dpi = 300, bbox_inches='tight')
adata_spatial.write_h5ad("/home/data/sdzl14/NSCLC/zong/spatial/adata_spatial.h5ad")
adata_spatial = sc.read_h5ad("/home/data/sdzl14/NSCLC/zong/spatial/adata_spatial.h5ad")
adata_spatial = adata_spatial.copy()
adata_spatial
ad = adata_spatial[adata_spatial.obs['library_id'] == 'ST_1.1.1']
sc.pl.spatial(ad,img_key="hires",color=['Cancer_P1/2/3'],edges=False,library_id='ST_1.1.1')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/tumor_ST_1.1.1.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/tumor_ST_1.1.1.png",dpi = 300, bbox_inches='tight')
ad = adata_spatial[adata_spatial.obs['library_id'] == 'ST_1.2.1']
sc.pl.spatial(ad,img_key="hires",color=['Neutrophils'],edges=False,library_id='ST_1.2.1')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/Macrophage_ST_1.2.1.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/Macrophage_ST_1.2.1.png",dpi = 300, bbox_inches='tight')
ad = adata_spatial[adata_spatial.obs['library_id'] == 'ST_1.1.1']
sc.pl.spatial(ad,img_key="hires",color=['Neutrophils'],edges=False,library_id='ST_1.1.1')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/Macrophage_ST_1.1.1.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/Macrophage_ST_1.1.1.png",dpi = 300, bbox_inches='tight')
ad = adata_spatial[adata_spatial.obs['library_id'] == 'ST_1.4.1']
sc.pl.spatial(ad,img_key="hires",color=['Neutrophils'],edges=False,library_id='ST_1.4.1')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/Macrophage_ST_1.4.1.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/Macrophage_ST_1.4.1.png",dpi = 300, bbox_inches='tight')
# 步骤2: 获取 CHI3L1 表达量（使用 log1p 标准化后的数据）
# 步骤1: 计算乘积
cancer_score1 = adata_spatial.obs['Cancer_P1/2/3']
cancer_score2 = adata_spatial.obs['Cancer_P4']
cancer_score3 = adata_spatial.obs['Cancer_P5']
cancer_score = cancer_score1 + cancer_score2 + cancer_score3
macrophage_score = adata_spatial.obs['Macrophage']
Fibroblast_score = adata_spatial.obs['Fibroblasts']
macrophage_cancer = cancer_score * macrophage_score

# 步骤2: 添加到 obs
adata_spatial.obs['Macrophage_Cancer_P1/2/3'] = macrophage_cancer

# 步骤1: 获取原始值
raw_values = adata_spatial.obs['Fibroblasts'].values

# 步骤2: 对数变换减小极端值影响（可选但推荐）
# 添加小常数避免log(0)错误
import numpy as np
log_transformed = np.log1p(raw_values)  # log(x+1)

# 步骤3: Min-Max标准化到[0,1]范围
min_val = np.min(log_transformed)
max_val = np.max(log_transformed)

# 避免除零错误
if max_val - min_val > 1e-10:
    scaled_values = (log_transformed - min_val) / (max_val - min_val)
else:
    scaled_values = np.zeros_like(log_transformed)  # 全零处理

# 步骤4: 添加标准化后的列
adata_spatial.obs['Fibroblasts'] = scaled_values


POSTN_expr = adata_spatial[:, 'POSTN'].X.toarray().flatten()
ADH1B_expr = adata_spatial[:, 'ADH1B'].X.toarray().flatten()
COL15A1_expr = adata_spatial[:, 'COL15A1'].X.toarray().flatten()
MLSN_expr = adata_spatial[:, 'MLSN'].X.toarray().flatten()
PI16_expr = adata_spatial[:, 'PI16'].X.toarray().flatten()
LAMA2_expr = adata_spatial[:, 'LAMA2'].X.toarray().flatten()
COL11A1_expr = adata_spatial[:, 'COL11A1'].X.toarray().flatten()

cancer_NME2 = NME2_expr * cancer_score
cancer_IGLC2  = IGLC2_expr * cancer_score
cancer_CLU = CLU_expr * cancer_score
cancer_ALDOA  = ALDOA_expr * cancer_score
Fibro_cancer_NME2 = cancer_NME2*macrophage_score
Fibro_cancer_IGLC2 = cancer_IGLC2*macrophage_score
Fibro_cancer_CLU = cancer_CLU*macrophage_score
Fibro_cancer_ALDOA = cancer_ALDOA*macrophage_score
adata_spatial.obs['Fibro_cancer_NME2'] = Fibro_cancer_NME2
adata_spatial.obs['Fibro_cancer_IGLC2'] = Fibro_cancer_IGLC2
adata_spatial.obs['Fibro_cancer_CLU'] = Fibro_cancer_CLU
adata_spatial.obs['Fibro_cancer_ALDOA'] = Fibro_cancer_ALDOA



# 步骤1: 获取原始值
raw_values = adata_spatial.obs['Fibro_cancer_NME2'].values

# 步骤2: 对数变换减小极端值影响（可选但推荐）
# 添加小常数避免log(0)错误
import numpy as np
log_transformed = np.log1p(raw_values)  # log(x+1)

# 步骤3: Min-Max标准化到[0,1]范围
min_val = np.min(log_transformed)
max_val = np.max(log_transformed)

# 避免除零错误
if max_val - min_val > 1e-10:
    scaled_values = (log_transformed - min_val) / (max_val - min_val)
else:
    scaled_values = np.zeros_like(log_transformed)  # 全零处理

# 步骤4: 添加标准化后的列
adata_spatial.obs['Fibro_cancer_NME2'] = scaled_values



    # 步骤4: 添加到 obs

ad = adata_spatial[adata_spatial.obs['library_id'] == 'ST_1.1.1']

sc.pl.spatial(
        ad,
        img_key="hires",
        color=['Macrophage_Cancer_P1/2/3'],
        edges=False,
        title="Macrophage_Cancer_P1/2/3 Score",
        cmap='coolwarm',library_id='ST_1.1.1')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/ST_1.1.1Macrophage_Cancer_.png", dpi=300)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/ST_1.1.1Macrophage_Cancer_.pdf", dpi=300)
ad = adata_spatial[adata_spatial.obs['library_id'] == 'ST_1.2.1']

sc.pl.spatial(
        ad,
        img_key="hires",
        color=['Macrophage_Cancer_P1/2/3'],
        edges=False,
        title="Macrophage_Cancer_P1/2/3 Score",
        cmap='coolwarm',library_id='ST_1.2.1')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/ST_1.2.1Macrophage_Cancer_.png", dpi=300)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/ST_1.2.1Macrophage_Cancer_.pdf", dpi=300)
ad = adata_spatial[adata_spatial.obs['library_id'] == 'ST_1.3.1']

sc.pl.spatial(
        ad,
        img_key="hires",
        color=['Macrophage_Cancer_P1/2/3'],
        edges=False,
        title="Macrophage_Cancer_P1/2/3 Score",
        cmap='coolwarm',library_id='ST_1.3.1')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/ST_1.3.1Macrophage_Cancer_.png", dpi=300)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/ST_1.3.1Macrophage_Cancer_.pdf", dpi=300)
lib_ids =['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1', 'ST_2.1.1',
       'ST_2.2.1', 'ST_2.2.2', 'ST_3.1.1', 'ST_3.2.1', 'ST_3.2.2', 'ST_3.3.1',
       'ST_5.2.1', 'ST_5.2.2', 'ST_5.3.1', 'ST_6.2.1', 'ST_6.3.1', 'ST_6.4.1']
lib_ids =['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2']

for lib_id in lib_ids:
    try:
        ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
        
        # 确保列存在
        if 'Macro_cancer_CLU' not in ad.obs:
            print(f"Warning: 'Macro_cancer_CLU' not found in {lib_id}")
            continue
            
        sc.pl.spatial(
            ad,
            img_key="hires",
            color=['Macro_cancer_CLU'],
            edges=False,
            title=f"Macro_cancer_CLU in {lib_id}",
            cmap='coolwarm',
            library_id=lib_id,
            show=False,
            vmin=0,  # 设置最小值
            vmax=2   # 设置最大值
        )
        
        safe_lib_id = lib_id.replace('.', '_')
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_Macro_cancer_CLU.png", dpi=300)
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_Macro_cancer_CLU.pdf", dpi=300)
        plt.close()
        
    except Exception as e:
        print(f"Error processing {lib_id}: {str(e)}")
        plt.close()  # 确保关闭可能残留的图形
# 使用squidpy创建更专业的可视化

sample_id = 'ST_1.1.1'
adata_sample = adata_spatial[adata_spatial.obs['library_id'] == sample_id].copy()

# 选择前几个主要细胞类型进行展示
main_cell_types = [ 'Macrophage', 'CD8+T', 'CD4+T', 'Mast cells','Monocytes','DC']



# 创建子图
fig, axes = plt.subplots(1, len(main_cell_types), figsize=(30, 4))

for idx, cell_type in enumerate(main_cell_types):
    sq.pl.spatial_scatter(adata_sample, 
                         color=f'{cell_type}', 
                         size=2, 
                         img=True, 
                         ax=axes[idx],
                         library_id=sample_id,
                         title=cell_type)
    
plt.tight_layout()
plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/main_cell_types_spatial_{sample_id}.pdf", 
            dpi=300, bbox_inches='tight')
import squidpy as sq
import matplotlib.pyplot as plt

# 选择一个样本进行展示（例如ST_1.1.1）
sample_id = 'ST_1.3.1'
adata_sample = adata_spatial[adata_spatial.obs['library_id'] == sample_id].copy()

# 选择前几个主要细胞类型进行展示
main_cell_types = ['Macrophage', 'CD8+T', 'CD4+T', 'Mast cells', 'Monocytes', 'DC']

# 获取所有细胞类型的表达量
expression_values = []
for cell_type in main_cell_types:
    expression_values.extend(adata_sample.obs[cell_type].values)

# 找到最大值和最小值
vmin = min(expression_values)
vmax = max(expression_values)

# 创建子图
fig, axes = plt.subplots(1, len(main_cell_types), figsize=(30, 4))

for idx, cell_type in enumerate(main_cell_types):
    sq.pl.spatial_scatter(
        adata_sample,
        color=f'{cell_type}',
        size=2,
        img=True,
        ax=axes[idx],
        library_id=sample_id,
        title=cell_type,
        vmin=vmin,  # 设置统一的最小值
        vmax=vmax   # 设置统一的最大值
    )

plt.tight_layout()
plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/main_cell_types_spatial_{sample_id}.pdf", 
            dpi=300, bbox_inches='tight')
plt.show()
# 定义缺氧适应基因列表
hypoxia_adaptation_list = [
    "VEGFA",     # 血管新生（angiogenesis）
    "SLC2A1",    # GLUT1，葡萄糖摄取
    "LDHA",      # 乳酸代谢
    "ENO1",      # 糖酵解
    "PGK1",      # 糖酵解
    "HK2",       # 糖酵解启动酶
    "ALDOA",     # 糖酵解
    "PKM",       # 糖酵解
    "PFKFB3",    # 糖酵解调控
    "CA9",       # pH稳态调节
    "BNIP3",     # 缺氧诱导的自噬与凋亡调节
    "ADM",       # Adrenomedullin，促进血管舒张与存活
    "ANGPTL4",   # 调控血管通透性和代谢
    "DDIT4",     # REDD1，mTOR抑制，能量节省
    "EGLN3",     # PHD3，负反馈调控HIF1稳定性
    "P4HA1",     # 胶原羟化酶，促进基质稳定性
    "PLOD2",     # 胶原交联
    "LOX",       # 胶原交联酶
    "CXCL12",    # 趋化因子，与肿瘤侵袭有关
    "MCT4",      # 乳酸转运（SLC16A3）
    "NDRG1",     # 应激/代谢调节
    "SERPINE1",  # PAI-1，血管稳态与肿瘤相关
]

# 检查哪些基因在数据中存在
available_genes = [gene for gene in hypoxia_adaptation_list if gene in adata_spatial.var_names]
print(f"在数据中找到的基因: {available_genes}")
print(f"缺失的基因: {set(hypoxia_adaptation_list) - set(available_genes)}")

# 提取这些基因的表达值
if len(available_genes) > 0:
    # 获取基因表达矩阵（使用原始counts或者已经标准化的数据）
    gene_expression = adata_spatial[:, available_genes].X
    
    # 转换为密集矩阵（如果需要）
    if not isinstance(gene_expression, np.ndarray):
        gene_expression = gene_expression.toarray()
    
    # 计算每个spot的基因表达总和
    hypoxia_score = np.sum(gene_expression, axis=1)
    
    # 标准化处理
    # 方法1: Z-score标准化（均值为0，标准差为1）
    hypoxia_score_zscore = (hypoxia_score - np.mean(hypoxia_score)) / np.std(hypoxia_score)
    
    # 方法2: Min-Max标准化（0-1范围）
    min_val = np.min(hypoxia_score)
    max_val = np.max(hypoxia_score)
    if max_val - min_val > 1e-10:
        hypoxia_score_minmax = (hypoxia_score - min_val) / (max_val - min_val)
    else:
        hypoxia_score_minmax = np.zeros_like(hypoxia_score)
    
    # 添加到obs中
    adata_spatial.obs['hypoxia_adaptation_score_zscore'] = hypoxia_score_zscore
    adata_spatial.obs['hypoxia_adaptation_score_minmax'] = hypoxia_score_minmax
    
    # 绘制每个样本的空间图
    library_ids = adata_spatial.obs['library_id'].unique()
    library_id = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
    # 使用Z-score标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_zscore',
                edges=False,
                title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
                cmap='viridis',  # 可以改为 'coolwarm', 'RdYlBu_r' 等其他colormap
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_hypoxia_adaptation_zscore.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_hypoxia_adaptation_zscore.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 使用Min-Max标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_minmax',
                edges=False,
                title=f'Hypoxia Adaptation Score (Min-Max) - {lib_id}',
                cmap='viridis',
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_hypoxia_adaptation_minmax.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_hypoxia_adaptation_minmax.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 也可以绘制UMAP图查看整体分布
    sc.pl.umap(adata_spatial, 
               color=['hypoxia_adaptation_score_zscore'], 
               edges=False, 
               cmap='viridis')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.png", 
               dpi=300, bbox_inches='tight')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.pdf", 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print("缺氧适应基因得分计算完成并已生成图表")
    print(f"Z-score范围: {np.min(hypoxia_score_zscore):.3f} 到 {np.max(hypoxia_score_zscore):.3f}")
    print(f"Min-Max范围: {np.min(hypoxia_score_minmax):.3f} 到 {np.max(hypoxia_score_minmax):.3f}")
    
else:
    print("在数据中未找到任何缺氧适应相关基因")
# 为Z-score得分设置统一颜色尺度
hypoxia_zscore_min = adata_spatial.obs['hypoxia_adaptation_score_zscore'].min()
hypoxia_zscore_max = adata_spatial.obs['hypoxia_adaptation_score_zscore'].max()

# 为Min-Max得分设置统一颜色尺度（这个应该是0-1，但为了保险起见还是计算一下）
hypoxia_minmax_min = adata_spatial.obs['hypoxia_adaptation_score_minmax'].min()
hypoxia_minmax_max = adata_spatial.obs['hypoxia_adaptation_score_minmax'].max()
lib_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
# 绘制统一颜色尺度的缺氧适应基因得分图
for lib_id in lib_ids:
    try:
        ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
        
        # Z-score标准化结果
        sc.pl.spatial(
            ad,
            img_key="hires",
            color='hypoxia_adaptation_score_zscore',
            edges=False,
            title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
            cmap='viridis',
            library_id=lib_id,
            show=False,
            vmin=hypoxia_zscore_min,
            vmax=hypoxia_zscore_max
        )
        
        safe_lib_id = lib_id.replace('.', '_')
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_hypoxia_adaptation_zscore_unified.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error processing hypoxia Z-score for {lib_id}: {str(e)}")
        plt.close()

# 如果需要，也可以为其他感兴趣的基因设置统一颜色尺度
genes_of_interest = ['CHI3L1', 'CCL18', 'CXCL3', 'CLU', 'ALDOA', 'NME2']
for gene in genes_of_interest:
    if gene in adata_spatial.var_names:
        # 计算全局最小值和最大值
        gene_global_min = adata_spatial[:, gene].X.min()
        gene_global_max = adata_spatial[:, gene].X.max()
        
        # 为每个library_id绘制统一颜色尺度的图
        for lib_id in lib_ids:
            try:
                ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
                
                sc.pl.spatial(
                    ad,
                    img_key="hires",
                    color=[gene],
                    edges=False,
                    title=f'{gene} in {lib_id}',
                    cmap='viridis',
                    library_id=lib_id,
                    show=False,
                    vmin=gene_global_min,
                    vmax=gene_global_max
                )
                
                safe_lib_id = lib_id.replace('.', '_')
                plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_{gene}_unified.png", 
                           dpi=300, bbox_inches='tight')
                plt.close()
                
            except Exception as e:
                print(f"Error processing {gene} for {lib_id}: {str(e)}")
                plt.close()
src_fak_signaling_gene_list = [
    # Core kinases
    "SRC",       # Src tyrosine kinase
    "PTK2",      # FAK (Focal Adhesion Kinase)

    # Adaptor and scaffold proteins
    "PXN",       # Paxillin
    "TLN1",      # Talin 1
    "VCL",       # Vinculin
    "ZYX",       # Zyxin
    "BCAR1",     # p130Cas
    "SHC1",      # SHC Adaptor Protein 1
    "GRB2",      # Growth factor receptor-bound protein 2
]
# 检查哪些基因在数据中存在
available_genes = [gene for gene in src_fak_signaling_gene_list if gene in adata_spatial.var_names]
print(f"在数据中找到的基因: {available_genes}")
print(f"缺失的基因: {set(src_fak_signaling_gene_list) - set(available_genes)}")

# 提取这些基因的表达值
if len(available_genes) > 0:
    # 获取基因表达矩阵（使用原始counts或者已经标准化的数据）
    gene_expression = adata_spatial[:, available_genes].X
    
    # 转换为密集矩阵（如果需要）
    if not isinstance(gene_expression, np.ndarray):
        gene_expression = gene_expression.toarray()
    
    # 计算每个spot的基因表达总和
    hypoxia_score = np.sum(gene_expression, axis=1)
    
    # 标准化处理
    # 方法1: Z-score标准化（均值为0，标准差为1）
    hypoxia_score_zscore = (hypoxia_score - np.mean(hypoxia_score)) / np.std(hypoxia_score)
    
    # 方法2: Min-Max标准化（0-1范围）
    min_val = np.min(hypoxia_score)
    max_val = np.max(hypoxia_score)
    if max_val - min_val > 1e-10:
        hypoxia_score_minmax = (hypoxia_score - min_val) / (max_val - min_val)
    else:
        hypoxia_score_minmax = np.zeros_like(hypoxia_score)
    
    # 添加到obs中
    adata_spatial.obs['hypoxia_adaptation_score_zscore'] = hypoxia_score_zscore
    adata_spatial.obs['hypoxia_adaptation_score_minmax'] = hypoxia_score_minmax
    
    # 绘制每个样本的空间图
    library_ids = adata_spatial.obs['library_id'].unique()
    library_id = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
    # 使用Z-score标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_zscore',
                edges=False,
                title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
                cmap='viridis',  # 可以改为 'coolwarm', 'RdYlBu_r' 等其他colormap
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_src_fak_signaling.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_src_fak_signaling.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 使用Min-Max标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_minmax',
                edges=False,
                title=f'Hypoxia Adaptation Score (Min-Max) - {lib_id}',
                cmap='viridis',
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_src_fak_signalingminmax.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_src_fak_signalingminmax.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 也可以绘制UMAP图查看整体分布
    sc.pl.umap(adata_spatial, 
               color=['hypoxia_adaptation_score_zscore'], 
               edges=False, 
               cmap='viridis')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.png", 
               dpi=300, bbox_inches='tight')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.pdf", 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print("缺氧适应基因得分计算完成并已生成图表")
    print(f"Z-score范围: {np.min(hypoxia_score_zscore):.3f} 到 {np.max(hypoxia_score_zscore):.3f}")
    print(f"Min-Max范围: {np.min(hypoxia_score_minmax):.3f} 到 {np.max(hypoxia_score_minmax):.3f}")
    
else:
    print("在数据中未找到任何缺氧适应相关基因")
# 为Z-score得分设置统一颜色尺度
hypoxia_zscore_min = adata_spatial.obs['hypoxia_adaptation_score_zscore'].min()
hypoxia_zscore_max = adata_spatial.obs['hypoxia_adaptation_score_zscore'].max()

# 为Min-Max得分设置统一颜色尺度（这个应该是0-1，但为了保险起见还是计算一下）
hypoxia_minmax_min = adata_spatial.obs['hypoxia_adaptation_score_minmax'].min()
hypoxia_minmax_max = adata_spatial.obs['hypoxia_adaptation_score_minmax'].max()
lib_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
# 绘制统一颜色尺度的缺氧适应基因得分图
for lib_id in lib_ids:
    try:
        ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
        
        # Z-score标准化结果
        sc.pl.spatial(
            ad,
            img_key="hires",
            color='hypoxia_adaptation_score_zscore',
            edges=False,
            title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
            cmap='viridis',
            library_id=lib_id,
            show=False,
            vmin=hypoxia_zscore_min,
            vmax=hypoxia_zscore_max
        )
        
        safe_lib_id = lib_id.replace('.', '_')
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_src_fak_signaling_unified.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error processing hypoxia Z-score for {lib_id}: {str(e)}")
        plt.close()

pseudopodia_cytoskeleton_gene_list = [
    # Rho family GTPases
    "RAC1", "CDC42", "RHOA",

    # WAVE complex
    "WASF2",  # WAVE2
    "NCKAP1",
    "ABI1",

    # ARP2/3 complex
    "ARPC1B",
    "ACTR2",
    "ACTR3"
]
# 检查哪些基因在数据中存在
available_genes = [gene for gene in pseudopodia_cytoskeleton_gene_list if gene in adata_spatial.var_names]
print(f"在数据中找到的基因: {available_genes}")
print(f"缺失的基因: {set(pseudopodia_cytoskeleton_gene_list) - set(available_genes)}")

# 提取这些基因的表达值
if len(available_genes) > 0:
    # 获取基因表达矩阵（使用原始counts或者已经标准化的数据）
    gene_expression = adata_spatial[:, available_genes].X
    
    # 转换为密集矩阵（如果需要）
    if not isinstance(gene_expression, np.ndarray):
        gene_expression = gene_expression.toarray()
    
    # 计算每个spot的基因表达总和
    hypoxia_score = np.sum(gene_expression, axis=1)
    
    # 标准化处理
    # 方法1: Z-score标准化（均值为0，标准差为1）
    hypoxia_score_zscore = (hypoxia_score - np.mean(hypoxia_score)) / np.std(hypoxia_score)
    
    # 方法2: Min-Max标准化（0-1范围）
    min_val = np.min(hypoxia_score)
    max_val = np.max(hypoxia_score)
    if max_val - min_val > 1e-10:
        hypoxia_score_minmax = (hypoxia_score - min_val) / (max_val - min_val)
    else:
        hypoxia_score_minmax = np.zeros_like(hypoxia_score)
    
    # 添加到obs中
    adata_spatial.obs['hypoxia_adaptation_score_zscore'] = hypoxia_score_zscore
    adata_spatial.obs['hypoxia_adaptation_score_minmax'] = hypoxia_score_minmax
    
    # 绘制每个样本的空间图
    library_ids = adata_spatial.obs['library_id'].unique()
    library_id = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
    # 使用Z-score标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_zscore',
                edges=False,
                title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
                cmap='viridis',  # 可以改为 'coolwarm', 'RdYlBu_r' 等其他colormap
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_pseudopodia_cytoskeleton.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_pseudopodia_cytoskeleton.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 使用Min-Max标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_minmax',
                edges=False,
                title=f'Hypoxia Adaptation Score (Min-Max) - {lib_id}',
                cmap='viridis',
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_pseudopodia_cytoskeletonminmax.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_pseudopodia_cytoskeletonminmax.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 也可以绘制UMAP图查看整体分布
    sc.pl.umap(adata_spatial, 
               color=['hypoxia_adaptation_score_zscore'], 
               edges=False, 
               cmap='viridis')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.png", 
               dpi=300, bbox_inches='tight')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.pdf", 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print("缺氧适应基因得分计算完成并已生成图表")
    print(f"Z-score范围: {np.min(hypoxia_score_zscore):.3f} 到 {np.max(hypoxia_score_zscore):.3f}")
    print(f"Min-Max范围: {np.min(hypoxia_score_minmax):.3f} 到 {np.max(hypoxia_score_minmax):.3f}")
    
else:
    print("在数据中未找到任何缺氧适应相关基因")
# 为Z-score得分设置统一颜色尺度
hypoxia_zscore_min = adata_spatial.obs['hypoxia_adaptation_score_zscore'].min()
hypoxia_zscore_max = adata_spatial.obs['hypoxia_adaptation_score_zscore'].max()

# 为Min-Max得分设置统一颜色尺度（这个应该是0-1，但为了保险起见还是计算一下）
hypoxia_minmax_min = adata_spatial.obs['hypoxia_adaptation_score_minmax'].min()
hypoxia_minmax_max = adata_spatial.obs['hypoxia_adaptation_score_minmax'].max()
lib_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
# 绘制统一颜色尺度的缺氧适应基因得分图
for lib_id in lib_ids:
    try:
        ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
        
        # Z-score标准化结果
        sc.pl.spatial(
            ad,
            img_key="hires",
            color='hypoxia_adaptation_score_zscore',
            edges=False,
            title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
            cmap='viridis',
            library_id=lib_id,
            show=False,
            vmin=hypoxia_zscore_min,
            vmax=hypoxia_zscore_max
        )
        
        safe_lib_id = lib_id.replace('.', '_')
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_pseudopodia_cytoskeletonminmax_unified.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error processing hypoxia Z-score for {lib_id}: {str(e)}")
        plt.close()
camp_pka_signaling_gene_list = [
    # cAMP synthesis and degradation
    "ADCY1", "ADCY2", "ADCY3", "ADCY4", "ADCY5", "ADCY6", "ADCY7", "ADCY8", "ADCY9",  # Adenylate cyclases
    "PDE4A", "PDE4B", "PDE4D", "PDE7A", "PDE7B",  # Phosphodiesterases (cAMP degradation)

    # cAMP effectors
    "PRKACA", "PRKACB", "PRKACG",    # PKA catalytic subunits
    "PRKAR1A", "PRKAR1B", "PRKAR2A", "PRKAR2B",  # PKA regulatory subunits

    # PKA downstream effectors (selected)
    "CREB1",     # cAMP response element-binding protein
    "ATF1",      # Activating transcription factor 1
    "BDNF",      # Brain-derived neurotrophic factor (PKA–CREB target)
    "NR4A1",     # Immediate early gene, CREB target
    "FOS", "JUN",  # AP-1 components, CREB-related
    "VASP",      # Actin regulatory protein, PKA substrate
    "PPP1R1A",   # Inhibitor-1, regulated by PKA

    # GPCR upstream regulators (optional, if interested in receptor-level input)
    "GNAS",      # Gαs subunit (activates adenylate cyclase)
    "GNAI1", "GNAI2", "GNAI3",  # Gαi subunits (inhibit AC)
    "GNB1", "GNG2"              # βγ subunits
]
available_genes = [gene for gene in camp_pka_signaling_gene_list if gene in adata_spatial.var_names]
print(f"在数据中找到的基因: {available_genes}")
print(f"缺失的基因: {set(camp_pka_signaling_gene_list) - set(available_genes)}")

# 提取这些基因的表达值
if len(available_genes) > 0:
    # 获取基因表达矩阵（使用原始counts或者已经标准化的数据）
    gene_expression = adata_spatial[:, available_genes].X
    
    # 转换为密集矩阵（如果需要）
    if not isinstance(gene_expression, np.ndarray):
        gene_expression = gene_expression.toarray()
    
    # 计算每个spot的基因表达总和
    hypoxia_score = np.sum(gene_expression, axis=1)
    
    # 标准化处理
    # 方法1: Z-score标准化（均值为0，标准差为1）
    hypoxia_score_zscore = (hypoxia_score - np.mean(hypoxia_score)) / np.std(hypoxia_score)
    
    # 方法2: Min-Max标准化（0-1范围）
    min_val = np.min(hypoxia_score)
    max_val = np.max(hypoxia_score)
    if max_val - min_val > 1e-10:
        hypoxia_score_minmax = (hypoxia_score - min_val) / (max_val - min_val)
    else:
        hypoxia_score_minmax = np.zeros_like(hypoxia_score)
    
    # 添加到obs中
    adata_spatial.obs['hypoxia_adaptation_score_zscore'] = hypoxia_score_zscore
    adata_spatial.obs['hypoxia_adaptation_score_minmax'] = hypoxia_score_minmax
    
    # 绘制每个样本的空间图
    library_ids = adata_spatial.obs['library_id'].unique()
    library_id = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
    # 使用Z-score标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_zscore',
                edges=False,
                title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
                cmap='viridis',  # 可以改为 'coolwarm', 'RdYlBu_r' 等其他colormap
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_camp_pka.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_camp_pka.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 使用Min-Max标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_minmax',
                edges=False,
                title=f'Hypoxia Adaptation Score (Min-Max) - {lib_id}',
                cmap='viridis',
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_camp_pkaminmax.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_camp_pkaminmax.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 也可以绘制UMAP图查看整体分布
    sc.pl.umap(adata_spatial, 
               color=['hypoxia_adaptation_score_zscore'], 
               edges=False, 
               cmap='viridis')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.png", 
               dpi=300, bbox_inches='tight')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.pdf", 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print("缺氧适应基因得分计算完成并已生成图表")
    print(f"Z-score范围: {np.min(hypoxia_score_zscore):.3f} 到 {np.max(hypoxia_score_zscore):.3f}")
    print(f"Min-Max范围: {np.min(hypoxia_score_minmax):.3f} 到 {np.max(hypoxia_score_minmax):.3f}")
    
else:
    print("在数据中未找到任何缺氧适应相关基因")
# 为Z-score得分设置统一颜色尺度
hypoxia_zscore_min = adata_spatial.obs['hypoxia_adaptation_score_zscore'].min()
hypoxia_zscore_max = adata_spatial.obs['hypoxia_adaptation_score_zscore'].max()

# 为Min-Max得分设置统一颜色尺度（这个应该是0-1，但为了保险起见还是计算一下）
hypoxia_minmax_min = adata_spatial.obs['hypoxia_adaptation_score_minmax'].min()
hypoxia_minmax_max = adata_spatial.obs['hypoxia_adaptation_score_minmax'].max()
lib_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
# 绘制统一颜色尺度的缺氧适应基因得分图
for lib_id in lib_ids:
    try:
        ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
        
        # Z-score标准化结果
        sc.pl.spatial(
            ad,
            img_key="hires",
            color='hypoxia_adaptation_score_zscore',
            edges=False,
            title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
            cmap='viridis',
            library_id=lib_id,
            show=False,
            vmin=hypoxia_zscore_min,
            vmax=hypoxia_zscore_max
        )
        
        safe_lib_id = lib_id.replace('.', '_')
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_camp_pkaunified.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error processing hypoxia Z-score for {lib_id}: {str(e)}")
        plt.close()
type_I_IFN_response_markers = [
    "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
    "IFI6", "IFI27", "IFIT1", "IFIT2", "IFIT3",
    "STAT1", "IRF7", "IRF9", 
    "CXCL10", "CXCL9"
]
available_genes = [gene for gene in type_I_IFN_response_markers if gene in adata_spatial.var_names]
print(f"在数据中找到的基因: {available_genes}")
print(f"缺失的基因: {set(type_I_IFN_response_markers) - set(available_genes)}")

# 提取这些基因的表达值
if len(available_genes) > 0:
    # 获取基因表达矩阵（使用原始counts或者已经标准化的数据）
    gene_expression = adata_spatial[:, available_genes].X
    
    # 转换为密集矩阵（如果需要）
    if not isinstance(gene_expression, np.ndarray):
        gene_expression = gene_expression.toarray()
    
    # 计算每个spot的基因表达总和
    hypoxia_score = np.sum(gene_expression, axis=1)
    
    # 标准化处理
    # 方法1: Z-score标准化（均值为0，标准差为1）
    hypoxia_score_zscore = (hypoxia_score - np.mean(hypoxia_score)) / np.std(hypoxia_score)
    
    # 方法2: Min-Max标准化（0-1范围）
    min_val = np.min(hypoxia_score)
    max_val = np.max(hypoxia_score)
    if max_val - min_val > 1e-10:
        hypoxia_score_minmax = (hypoxia_score - min_val) / (max_val - min_val)
    else:
        hypoxia_score_minmax = np.zeros_like(hypoxia_score)
    
    # 添加到obs中
    adata_spatial.obs['hypoxia_adaptation_score_zscore'] = hypoxia_score_zscore
    adata_spatial.obs['hypoxia_adaptation_score_minmax'] = hypoxia_score_minmax
    
    # 绘制每个样本的空间图
    library_ids = adata_spatial.obs['library_id'].unique()
    library_id = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
    # 使用Z-score标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_zscore',
                edges=False,
                title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
                cmap='viridis',  # 可以改为 'coolwarm', 'RdYlBu_r' 等其他colormap
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_type_I_IFN.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_type_I_IFN.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 使用Min-Max标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_minmax',
                edges=False,
                title=f'Hypoxia Adaptation Score (Min-Max) - {lib_id}',
                cmap='viridis',
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_type_I_IFNminmax.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_type_I_IFNminmax.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 也可以绘制UMAP图查看整体分布
    sc.pl.umap(adata_spatial, 
               color=['hypoxia_adaptation_score_zscore'], 
               edges=False, 
               cmap='viridis')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.png", 
               dpi=300, bbox_inches='tight')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.pdf", 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print("缺氧适应基因得分计算完成并已生成图表")
    print(f"Z-score范围: {np.min(hypoxia_score_zscore):.3f} 到 {np.max(hypoxia_score_zscore):.3f}")
    print(f"Min-Max范围: {np.min(hypoxia_score_minmax):.3f} 到 {np.max(hypoxia_score_minmax):.3f}")
    
else:
    print("在数据中未找到任何缺氧适应相关基因")
# 为Z-score得分设置统一颜色尺度
hypoxia_zscore_min = adata_spatial.obs['hypoxia_adaptation_score_zscore'].min()
hypoxia_zscore_max = adata_spatial.obs['hypoxia_adaptation_score_zscore'].max()

# 为Min-Max得分设置统一颜色尺度（这个应该是0-1，但为了保险起见还是计算一下）
hypoxia_minmax_min = adata_spatial.obs['hypoxia_adaptation_score_minmax'].min()
hypoxia_minmax_max = adata_spatial.obs['hypoxia_adaptation_score_minmax'].max()
lib_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
# 绘制统一颜色尺度的缺氧适应基因得分图
for lib_id in lib_ids:
    try:
        ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
        
        # Z-score标准化结果
        sc.pl.spatial(
            ad,
            img_key="hires",
            color='hypoxia_adaptation_score_zscore',
            edges=False,
            title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
            cmap='viridis',
            library_id=lib_id,
            show=False,
            vmin=hypoxia_zscore_min,
            vmax=hypoxia_zscore_max
        )
        
        safe_lib_id = lib_id.replace('.', '_')
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_type_I_IFNunified.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error processing hypoxia Z-score for {lib_id}: {str(e)}")
        plt.close()
PI3K_AKT_mTOR_genes = [
    
    # PI3K 复合体
    "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG",  # Catalytic subunits
    "PIK3R1", "PIK3R2", "PIK3R3",            # Regulatory subunits

    # PTEN 抑制因子
    "PTEN",  # tumor suppressor

    # AKT 家族
    "AKT1", "AKT2", "AKT3",

    # mTOR 路径关键因子
    "MTOR", "RICTOR", "RPTOR",  # mTORC1 & mTORC2 components
    "MLST8", "DEPTOR",

]
available_genes = [gene for gene in PI3K_AKT_mTOR_genes if gene in adata_spatial.var_names]
print(f"在数据中找到的基因: {available_genes}")
print(f"缺失的基因: {set(PI3K_AKT_mTOR_genes) - set(available_genes)}")

# 提取这些基因的表达值
if len(available_genes) > 0:
    # 获取基因表达矩阵（使用原始counts或者已经标准化的数据）
    gene_expression = adata_spatial[:, available_genes].X
    
    # 转换为密集矩阵（如果需要）
    if not isinstance(gene_expression, np.ndarray):
        gene_expression = gene_expression.toarray()
    
    # 计算每个spot的基因表达总和
    hypoxia_score = np.sum(gene_expression, axis=1)
    
    # 标准化处理
    # 方法1: Z-score标准化（均值为0，标准差为1）
    hypoxia_score_zscore = (hypoxia_score - np.mean(hypoxia_score)) / np.std(hypoxia_score)
    
    # 方法2: Min-Max标准化（0-1范围）
    min_val = np.min(hypoxia_score)
    max_val = np.max(hypoxia_score)
    if max_val - min_val > 1e-10:
        hypoxia_score_minmax = (hypoxia_score - min_val) / (max_val - min_val)
    else:
        hypoxia_score_minmax = np.zeros_like(hypoxia_score)
    
    # 添加到obs中
    adata_spatial.obs['hypoxia_adaptation_score_zscore'] = hypoxia_score_zscore
    adata_spatial.obs['hypoxia_adaptation_score_minmax'] = hypoxia_score_minmax
    
    # 绘制每个样本的空间图
    library_ids = adata_spatial.obs['library_id'].unique()
    library_id = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
    # 使用Z-score标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_zscore',
                edges=False,
                title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
                cmap='viridis',  # 可以改为 'coolwarm', 'RdYlBu_r' 等其他colormap
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_PI3K_AKT.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_PI3K_AKT.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 使用Min-Max标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_minmax',
                edges=False,
                title=f'Hypoxia Adaptation Score (Min-Max) - {lib_id}',
                cmap='viridis',
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_PI3K_AKTminmax.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_PI3K_AKTminmax.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 也可以绘制UMAP图查看整体分布
    sc.pl.umap(adata_spatial, 
               color=['hypoxia_adaptation_score_zscore'], 
               edges=False, 
               cmap='viridis')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.png", 
               dpi=300, bbox_inches='tight')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.pdf", 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print("缺氧适应基因得分计算完成并已生成图表")
    print(f"Z-score范围: {np.min(hypoxia_score_zscore):.3f} 到 {np.max(hypoxia_score_zscore):.3f}")
    print(f"Min-Max范围: {np.min(hypoxia_score_minmax):.3f} 到 {np.max(hypoxia_score_minmax):.3f}")
    
else:
    print("在数据中未找到任何缺氧适应相关基因")
# 为Z-score得分设置统一颜色尺度
hypoxia_zscore_min = adata_spatial.obs['hypoxia_adaptation_score_zscore'].min()
hypoxia_zscore_max = adata_spatial.obs['hypoxia_adaptation_score_zscore'].max()

# 为Min-Max得分设置统一颜色尺度（这个应该是0-1，但为了保险起见还是计算一下）
hypoxia_minmax_min = adata_spatial.obs['hypoxia_adaptation_score_minmax'].min()
hypoxia_minmax_max = adata_spatial.obs['hypoxia_adaptation_score_minmax'].max()
lib_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
# 绘制统一颜色尺度的缺氧适应基因得分图
for lib_id in lib_ids:
    try:
        ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
        
        # Z-score标准化结果
        sc.pl.spatial(
            ad,
            img_key="hires",
            color='hypoxia_adaptation_score_zscore',
            edges=False,
            title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
            cmap='viridis',
            library_id=lib_id,
            show=False,
            vmin=hypoxia_zscore_min,
            vmax=hypoxia_zscore_max
        )
        
        safe_lib_id = lib_id.replace('.', '_')
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}__PI3K_AKTunified.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error processing hypoxia Z-score for {lib_id}: {str(e)}")
        plt.close()
NFkB_pathway_genes = [
    # 受体（可引发NF-κB活化）
    "TNFRSF1A", "TNFRSF1B",  # TNF receptors
    "IL1R1", "IL1RAP",       # IL-1 receptors
    "TLR2", "TLR4", "TLR9",  # Toll-like receptors
    "CD40", "CD27", "CD30",

    # IKK复合物成员（激酶激活复合物）
    "CHUK",    # IKKα
    "IKBKB",   # IKKβ
    "IKBKG",   # NEMO/IKKγ

    # IκB 蛋白（NF-κB 抑制蛋白）
    "NFKBIA",  # IκBα
    "NFKBIB", "NFKBIE", "BCL3",

    # 核内转录因子（NF-κB 家族核心）
    "RELA",    # p65
    "NFKB1",   # p50
    "NFKB2",   # p52
    "RELB",    # RelB
    "REL",     # c-Rel

]
available_genes = [gene for gene in NFkB_pathway_genes if gene in adata_spatial.var_names]
print(f"在数据中找到的基因: {available_genes}")
print(f"缺失的基因: {set(NFkB_pathway_genes) - set(available_genes)}")

# 提取这些基因的表达值
if len(available_genes) > 0:
    # 获取基因表达矩阵（使用原始counts或者已经标准化的数据）
    gene_expression = adata_spatial[:, available_genes].X
    
    # 转换为密集矩阵（如果需要）
    if not isinstance(gene_expression, np.ndarray):
        gene_expression = gene_expression.toarray()
    
    # 计算每个spot的基因表达总和
    hypoxia_score = np.sum(gene_expression, axis=1)
    
    # 标准化处理
    # 方法1: Z-score标准化（均值为0，标准差为1）
    hypoxia_score_zscore = (hypoxia_score - np.mean(hypoxia_score)) / np.std(hypoxia_score)
    
    # 方法2: Min-Max标准化（0-1范围）
    min_val = np.min(hypoxia_score)
    max_val = np.max(hypoxia_score)
    if max_val - min_val > 1e-10:
        hypoxia_score_minmax = (hypoxia_score - min_val) / (max_val - min_val)
    else:
        hypoxia_score_minmax = np.zeros_like(hypoxia_score)
    
    # 添加到obs中
    adata_spatial.obs['hypoxia_adaptation_score_zscore'] = hypoxia_score_zscore
    adata_spatial.obs['hypoxia_adaptation_score_minmax'] = hypoxia_score_minmax
    
    # 绘制每个样本的空间图
    library_ids = adata_spatial.obs['library_id'].unique()
    library_id = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
    # 使用Z-score标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_zscore',
                edges=False,
                title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
                cmap='viridis',  # 可以改为 'coolwarm', 'RdYlBu_r' 等其他colormap
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_NFkB.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_NFkB.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 使用Min-Max标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_minmax',
                edges=False,
                title=f'Hypoxia Adaptation Score (Min-Max) - {lib_id}',
                cmap='viridis',
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_NFkBminmax.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_NFkBminmax.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 也可以绘制UMAP图查看整体分布
    sc.pl.umap(adata_spatial, 
               color=['hypoxia_adaptation_score_zscore'], 
               edges=False, 
               cmap='viridis')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.png", 
               dpi=300, bbox_inches='tight')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.pdf", 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print("缺氧适应基因得分计算完成并已生成图表")
    print(f"Z-score范围: {np.min(hypoxia_score_zscore):.3f} 到 {np.max(hypoxia_score_zscore):.3f}")
    print(f"Min-Max范围: {np.min(hypoxia_score_minmax):.3f} 到 {np.max(hypoxia_score_minmax):.3f}")
    
else:
    print("在数据中未找到任何缺氧适应相关基因")
# 为Z-score得分设置统一颜色尺度
hypoxia_zscore_min = adata_spatial.obs['hypoxia_adaptation_score_zscore'].min()
hypoxia_zscore_max = adata_spatial.obs['hypoxia_adaptation_score_zscore'].max()

# 为Min-Max得分设置统一颜色尺度（这个应该是0-1，但为了保险起见还是计算一下）
hypoxia_minmax_min = adata_spatial.obs['hypoxia_adaptation_score_minmax'].min()
hypoxia_minmax_max = adata_spatial.obs['hypoxia_adaptation_score_minmax'].max()
lib_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
# 绘制统一颜色尺度的缺氧适应基因得分图
for lib_id in lib_ids:
    try:
        ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
        
        # Z-score标准化结果
        sc.pl.spatial(
            ad,
            img_key="hires",
            color='hypoxia_adaptation_score_zscore',
            edges=False,
            title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
            cmap='viridis',
            library_id=lib_id,
            show=False,
            vmin=hypoxia_zscore_min,
            vmax=hypoxia_zscore_max
        )
        
        safe_lib_id = lib_id.replace('.', '_')
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_NFkBunified.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error processing hypoxia Z-score for {lib_id}: {str(e)}")
        plt.close()
e2f_genes = [
    # 核心E2F转录因子
    "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8",

    # E2F靶基因，参与DNA复制和细胞周期
    "CCNE1", "CCNE2",         # Cyclin E1/2 (E2F targets)
    "CDK2",                  # Cyclin-dependent kinase
    "CDC6", "CDC25A",        # DNA replication licensing
    "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",  # Mini-chromosome maintenance complex
    "PCNA",                  # Proliferation marker
    "TK1",                   # Thymidine kinase 1
    "RRM1", "RRM2",          # Ribonucleotide reductase subunits
    "DHFR",                  # Dihydrofolate reductase
    "CDT1",                  # DNA replication licensing factor
    "GMNN",                  # Geminin
    "ORC1", "ORC6",          # Origin recognition complex
    "CHEK1", "CHEK2",        # Checkpoint kinases (sometimes E2F-regulated)
    "BRCA1", "RAD51"         # DNA repair & S phase progression
]
available_genes = [gene for gene in e2f_genes if gene in adata_spatial.var_names]
print(f"在数据中找到的基因: {available_genes}")
print(f"缺失的基因: {set(e2f_genes) - set(available_genes)}")

# 提取这些基因的表达值
if len(available_genes) > 0:
    # 获取基因表达矩阵（使用原始counts或者已经标准化的数据）
    gene_expression = adata_spatial[:, available_genes].X
    
    # 转换为密集矩阵（如果需要）
    if not isinstance(gene_expression, np.ndarray):
        gene_expression = gene_expression.toarray()
    
    # 计算每个spot的基因表达总和
    hypoxia_score = np.sum(gene_expression, axis=1)
    
    # 标准化处理
    # 方法1: Z-score标准化（均值为0，标准差为1）
    hypoxia_score_zscore = (hypoxia_score - np.mean(hypoxia_score)) / np.std(hypoxia_score)
    
    # 方法2: Min-Max标准化（0-1范围）
    min_val = np.min(hypoxia_score)
    max_val = np.max(hypoxia_score)
    if max_val - min_val > 1e-10:
        hypoxia_score_minmax = (hypoxia_score - min_val) / (max_val - min_val)
    else:
        hypoxia_score_minmax = np.zeros_like(hypoxia_score)
    
    # 添加到obs中
    adata_spatial.obs['hypoxia_adaptation_score_zscore'] = hypoxia_score_zscore
    adata_spatial.obs['hypoxia_adaptation_score_minmax'] = hypoxia_score_minmax
    
    # 绘制每个样本的空间图
    library_ids = adata_spatial.obs['library_id'].unique()
    library_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
    # 使用Z-score标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_zscore',
                edges=False,
                title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
                cmap='viridis',  # 可以改为 'coolwarm', 'RdYlBu_r' 等其他colormap
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_e2f_genes.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_e2f_genes.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 使用Min-Max标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_minmax',
                edges=False,
                title=f'e2f_genes Score (Min-Max) - {lib_id}',
                cmap='viridis',
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_e2f_genesminmax.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_e2f_genesminmax.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    
    print("缺氧适应基因得分计算完成并已生成图表")
    print(f"Z-score范围: {np.min(hypoxia_score_zscore):.3f} 到 {np.max(hypoxia_score_zscore):.3f}")
    print(f"Min-Max范围: {np.min(hypoxia_score_minmax):.3f} 到 {np.max(hypoxia_score_minmax):.3f}")
    
else:
    print("在数据中未找到任何缺氧适应相关基因")
# 为Z-score得分设置统一颜色尺度
hypoxia_zscore_min = adata_spatial.obs['hypoxia_adaptation_score_zscore'].min()
hypoxia_zscore_max = adata_spatial.obs['hypoxia_adaptation_score_zscore'].max()

# 为Min-Max得分设置统一颜色尺度（这个应该是0-1，但为了保险起见还是计算一下）
hypoxia_minmax_min = adata_spatial.obs['hypoxia_adaptation_score_minmax'].min()
hypoxia_minmax_max = adata_spatial.obs['hypoxia_adaptation_score_minmax'].max()
lib_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
# 绘制统一颜色尺度的缺氧适应基因得分图
for lib_id in lib_ids:
    try:
        ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
        
        # Z-score标准化结果
        sc.pl.spatial(
            ad,
            img_key="hires",
            color='hypoxia_adaptation_score_zscore',
            edges=False,
            title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
            cmap='viridis',
            library_id=lib_id,
            show=False,
            vmin=hypoxia_zscore_min,
            vmax=hypoxia_zscore_max
        )
        
        safe_lib_id = lib_id.replace('.', '_')
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_e2f_genesunified.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error processing hypoxia Z-score for {lib_id}: {str(e)}")
        plt.close()
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pysal.explore import esda
from libpysal.weights import lat2W
import seaborn as sns
from statsmodels.formula.api import ols
import statsmodels.api as sm
sc.tl.score_genes(adata_spatial, gene_list=cadherin_binding_up_genes, score_name='E2F_score')
from libpysal.weights import KNN
from libpysal.cg import KDTree

coords = adata_spatial.obsm['spatial']
knn = KNN.from_array(coords, k=6)  # k=6 为常用设置，可调整
import esda
adata_spatial.obs['tumor_score'] = adata_spatial.obs['Cancer_P1/2/3']  # 或 P2 / P3 视具体分析

model = ols('E2F_score ~ tumor_score', data=adata_spatial.obs).fit()
print(model.summary())
sns.lmplot(x='tumor_score', y='E2F_score', data=adata_spatial.obs, scatter_kws={"s": 10})
plt.title('Spatial regression: Tumor vs E2F score')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/tumor_E2F_regression.png", 
          dpi=300, bbox_inches='tight')
# 按照 Tissue 分组计算回归并可视化
tissue_groups = adata_spatial.obs['Tissue'].unique()
n_tissues = len(tissue_groups)

# 创建一个字典来存储每个组织的回归结果
regression_results = {}

# 创建一个用于存储所有数据的 DataFrame
plot_data = adata_spatial.obs[['tumor_score', 'E2F_score', 'Tissue']].copy()

# 创建子图
fig, axes = plt.subplots(2, 2, figsize=(15, 12))
axes = axes.flatten()

# 为每个组织类型计算回归并绘制
for i, tissue in enumerate(tissue_groups):
    if i >= 4:  # 限制为4个子图
        break
        
    # 筛选当前组织类型的数据
    tissue_data = adata_spatial.obs[adata_spatial.obs['Tissue'] == tissue]
    
    # 拟合回归模型
    model = ols('E2F_score ~ tumor_score', data=tissue_data).fit()
    regression_results[tissue] = model
    
    # 打印回归摘要
    print(f"\n{Tissue} Regression Results:")
    print(model.summary())
    
    # 绘制散点图和回归线
    ax = axes[i]
    sns.scatterplot(data=tissue_data, x='tumor_score', y='E2F_score', s=10, alpha=0.7, ax=ax)
    
    # 添加回归线
    x_vals = np.array(ax.get_xlim())
    y_vals = model.params['Intercept'] + model.params['tumor_score'] * x_vals
    ax.plot(x_vals, y_vals, '--', color='red', linewidth=2)
    
    ax.set_title(f'{tissue} (R² = {model.rsquared:.3f})')
    ax.set_xlabel('Tumor Score')
    ax.set_ylabel('E2F Score')

# 删除多余的子图（如果有）
for i in range(len(tissue_groups), 4):
    fig.delaxes(axes[i])

plt.tight_layout()
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/tumor_cadherin_regression_by_tissue.png", 
            dpi=300, bbox_inches='tight')
plt.show()

# 创建一个组合图，显示所有组织类型的回归线比较
plt.figure(figsize=(10, 8))

# 使用不同颜色绘制所有组织类型的数据点和回归线
palette = sns.color_palette("husl", n_tissues)
for i, tissue in enumerate(tissue_groups):
    tissue_data = adata_spatial.obs[adata_spatial.obs['Tissue'] == tissue]
    
    # 绘制散点图
    sns.scatterplot(data=tissue_data, x='tumor_score', y='E2F_score', 
                    label=tissue, color=palette[i], s=10, alpha=0.7)
    
    # 拟合并绘制回归线
    model = ols('E2F_score ~ tumor_score', data=tissue_data).fit()
    x_vals = np.array(plt.xlim())
    y_vals = model.params['Intercept'] + model.params['tumor_score'] * x_vals
    plt.plot(x_vals, y_vals, '--', color=palette[i], linewidth=2, 
             label=f"{tissue} (R²={model.rsquared:.3f})")

plt.xlabel('Tumor Score')
plt.ylabel('E2F Score')
plt.title('E2F Score vs Tumor Score by Tissue Type')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/tumor_cadherin_regression_combined.png", 
            dpi=300, bbox_inches='tight')
plt.show()

# 创建回归统计摘要表
regression_summary = []
for tissue in tissue_groups:
    tissue_data = adata_spatial.obs[adata_spatial.obs['Tissue'] == tissue]
    model = ols('E2F_score ~ tumor_score', data=tissue_data).fit()
    
    regression_summary.append({
        'Tissue': tissue,
        'R_squared': model.rsquared,
        'R_squared_adj': model.rsquared_adj,
        'Coefficient': model.params['tumor_score'],
        'P_value': model.pvalues['tumor_score'],
        'Intercept': model.params['Intercept']
    })

regression_df = pd.DataFrame(regression_summary)
print("\nRegression Summary by Tissue:")
print(regression_df)

# 可视化回归统计量
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# R² 值比较
sns.barplot(data=regression_df, x='Tissue', y='R_squared', ax=axes[0])
axes[0].set_title('R² by Tissue Type')
axes[0].set_ylabel('R²')
axes[0].tick_params(axis='x', rotation=45)

# 回归系数比较
sns.barplot(data=regression_df, x='Tissue', y='Coefficient', ax=axes[1])
axes[1].set_title('Coefficient by Tissue Type')
axes[1].set_ylabel('Coefficient')
axes[1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/tumor_cadherin_regression_stats.png", 
            dpi=300, bbox_inches='tight')
# 创建带回归线的散点图
plt.figure(figsize=(12, 8))

# 使用seaborn的scatterplot绘制散点图
sns.scatterplot(data=adata_spatial.obs, 
                x='tumor_score', 
                y='E2F_score', 
                hue='Tissue', 
                s=5, 
                alpha=0.7)

# 为每种Tissue类型添加拟合线
tissue_groups = adata_spatial.obs['Tissue'].unique()
palette = sns.color_palette("husl", len(tissue_groups))

for i, tissue in enumerate(tissue_groups):
    tissue_data = adata_spatial.obs[adata_spatial.obs['Tissue'] == tissue]
    
    # 计算线性回归线
    slope, intercept = np.polyfit(tissue_data['tumor_score'], tissue_data['E2F_score'], 1)
    x_vals = np.array(plt.xlim())
    y_vals = intercept + slope * x_vals
    
    # 绘制回归线（使用较浅的颜色以区分原始数据点）
    plt.plot(x_vals, y_vals, '--', color=palette[i], linewidth=5, alpha=0.8)

plt.xlabel('Tumor Score')
plt.ylabel('E2F Score')
plt.title('E2F Score vs Tumor Score by Tissue Type with Regression Lines')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# 保存图像
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/tumor_cadherin_scatter_regression_by_tissue.png", 
            dpi=300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/tumor_cadherin_scatter_regression_by_tissue.pdf", 
            dpi=300, bbox_inches='tight')
hedgehog_signaling_genes = [
    # Hedgehog配体
    "SHH",     # Sonic Hedgehog
    "IHH",     # Indian Hedgehog
    "DHH",     # Desert Hedgehog

    # Hedgehog受体与调控因子
    "PTCH1",   # Patched 1 (主受体)
    "PTCH2",
    "SMO",     # Smoothened (信号转导关键蛋白)

    # 下游转录调控
    "GLI1", "GLI2", "GLI3",  # GLI转录因子家族
    "SUFU",     # 抑制因子 Suppressor of Fused

    # 正调控或靶基因
    "HHIP",     # Hedgehog-interacting protein, 负反馈靶基因
    "BCL2",     # 抗凋亡调控（GLI靶基因之一）
    "CCND1", "CCNE1",  # 细胞周期推进
    "FOXM1",    # 细胞增殖
    "JAG2",     # 非经典靶点，Notch互作

    # 与Hh信号协同的共因子
    "KIF7", "STK36", "GSK3B", "CSNK1G1",  # 信号复合物调控

    # 信号调控与反馈
    "ZIC2", "ZIC1",   # 调控GLI活性
    "SPOP",          # 负调控GLI降解
]
available_genes = [gene for gene in hedgehog_signaling_genes if gene in adata_spatial.var_names]
print(f"在数据中找到的基因: {available_genes}")
print(f"缺失的基因: {set(hedgehog_signaling_genes) - set(available_genes)}")

# 提取这些基因的表达值
if len(available_genes) > 0:
    # 获取基因表达矩阵（使用原始counts或者已经标准化的数据）
    gene_expression = adata_spatial[:, available_genes].X
    
    # 转换为密集矩阵（如果需要）
    if not isinstance(gene_expression, np.ndarray):
        gene_expression = gene_expression.toarray()
    
    # 计算每个spot的基因表达总和
    hypoxia_score = np.sum(gene_expression, axis=1)
    
    # 标准化处理
    # 方法1: Z-score标准化（均值为0，标准差为1）
    hypoxia_score_zscore = (hypoxia_score - np.mean(hypoxia_score)) / np.std(hypoxia_score)
    
    # 方法2: Min-Max标准化（0-1范围）
    min_val = np.min(hypoxia_score)
    max_val = np.max(hypoxia_score)
    if max_val - min_val > 1e-10:
        hypoxia_score_minmax = (hypoxia_score - min_val) / (max_val - min_val)
    else:
        hypoxia_score_minmax = np.zeros_like(hypoxia_score)
    
    # 添加到obs中
    adata_spatial.obs['hypoxia_adaptation_score_zscore'] = hypoxia_score_zscore
    adata_spatial.obs['hypoxia_adaptation_score_minmax'] = hypoxia_score_minmax
    
    # 绘制每个样本的空间图
    library_ids = adata_spatial.obs['library_id'].unique()
    library_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
    # 使用Z-score标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_zscore',
                edges=False,
                title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
                cmap='viridis',  # 可以改为 'coolwarm', 'RdYlBu_r' 等其他colormap
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_hedgehog.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_hedgehog.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 使用Min-Max标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_minmax',
                edges=False,
                title=f'Hypoxia Adaptation Score (Min-Max) - {lib_id}',
                cmap='viridis',
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_hedgehog_minmax.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_hedgehog_minmax.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 也可以绘制UMAP图查看整体分布
    sc.pl.umap(adata_spatial, 
               color=['hypoxia_adaptation_score_zscore'], 
               edges=False, 
               cmap='viridis')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.png", 
               dpi=300, bbox_inches='tight')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.pdf", 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print("缺氧适应基因得分计算完成并已生成图表")
    print(f"Z-score范围: {np.min(hypoxia_score_zscore):.3f} 到 {np.max(hypoxia_score_zscore):.3f}")
    print(f"Min-Max范围: {np.min(hypoxia_score_minmax):.3f} 到 {np.max(hypoxia_score_minmax):.3f}")
    
else:
    print("在数据中未找到任何缺氧适应相关基因")
# 为Z-score得分设置统一颜色尺度
hypoxia_zscore_min = adata_spatial.obs['hypoxia_adaptation_score_zscore'].min()
hypoxia_zscore_max = adata_spatial.obs['hypoxia_adaptation_score_zscore'].max()

# 为Min-Max得分设置统一颜色尺度（这个应该是0-1，但为了保险起见还是计算一下）
hypoxia_minmax_min = adata_spatial.obs['hypoxia_adaptation_score_minmax'].min()
hypoxia_minmax_max = adata_spatial.obs['hypoxia_adaptation_score_minmax'].max()
lib_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
# 绘制统一颜色尺度的缺氧适应基因得分图
for lib_id in lib_ids:
    try:
        ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
        
        # Z-score标准化结果
        sc.pl.spatial(
            ad,
            img_key="hires",
            color='hypoxia_adaptation_score_zscore',
            edges=False,
            title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
            cmap='viridis',
            library_id=lib_id,
            show=False,
            vmin=hypoxia_zscore_min,
            vmax=hypoxia_zscore_max
        )
        
        safe_lib_id = lib_id.replace('.', '_')
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_E2F2_unified.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error processing hypoxia Z-score for {lib_id}: {str(e)}")
        plt.close()
cadherin_binding_up_genes = [
    # EMT过程中典型上调的cadherins或其结合蛋白
    "CDH2",     # N-cadherin（取代CDH1，标志EMT）
    "CDH11",    # OB-cadherin，CAF常高表达
    "CTNNB1",   # β-catenin（WNT/EMT）
    "CTNND1",   # p120-catenin，调控黏附动态
    "JUP",      # plakoglobin（可增强细胞黏附）

    # EMT转录调控因子（调控 cadherin 表达）
    "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2",

    # cytoskeleton-linking and force-transducing
    "VCL",      # vinculin
    "ACTN4",    # α-actinin-4，促进细胞迁移
    "FERMT2",   # kindlin-2，CAF与invasion相关

    # Rho GTPases 及其调节因子（介导cadherin结合动力学）
    "RHOA", "RAC1", "CDC42",
    "ARHGEF7", "TIAM1",

    # 与CAF功能相关的integrin和ECM交互蛋白（间接影响cadherin binding）
    "ITGB1", "ITGA5",
    "LAMC1", "LAMB1",

    # desmosomal cadherins（上皮–间充质混合状态）
    "DSG2", "DSC2"
]
available_genes = [gene for gene in cadherin_binding_up_genes if gene in adata_spatial.var_names]
print(f"在数据中找到的基因: {available_genes}")
print(f"缺失的基因: {set(cadherin_binding_up_genes) - set(available_genes)}")

# 提取这些基因的表达值
if len(available_genes) > 0:
    # 获取基因表达矩阵（使用原始counts或者已经标准化的数据）
    gene_expression = adata_spatial[:, available_genes].X
    
    # 转换为密集矩阵（如果需要）
    if not isinstance(gene_expression, np.ndarray):
        gene_expression = gene_expression.toarray()
    
    # 计算每个spot的基因表达总和
    hypoxia_score = np.sum(gene_expression, axis=1)
    
    # 标准化处理
    # 方法1: Z-score标准化（均值为0，标准差为1）
    hypoxia_score_zscore = (hypoxia_score - np.mean(hypoxia_score)) / np.std(hypoxia_score)
    
    # 方法2: Min-Max标准化（0-1范围）
    min_val = np.min(hypoxia_score)
    max_val = np.max(hypoxia_score)
    if max_val - min_val > 1e-10:
        hypoxia_score_minmax = (hypoxia_score - min_val) / (max_val - min_val)
    else:
        hypoxia_score_minmax = np.zeros_like(hypoxia_score)
    
    # 添加到obs中
    adata_spatial.obs['hypoxia_adaptation_score_zscore'] = hypoxia_score_zscore
    adata_spatial.obs['hypoxia_adaptation_score_minmax'] = hypoxia_score_minmax
    
    # 绘制每个样本的空间图
    library_ids = adata_spatial.obs['library_id'].unique()
    library_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
    # 使用Z-score标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_zscore',
                edges=False,
                title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
                cmap='viridis',  # 可以改为 'coolwarm', 'RdYlBu_r' 等其他colormap
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_cadherin.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_cadherin.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 使用Min-Max标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_minmax',
                edges=False,
                title=f'Hypoxia Adaptation Score (Min-Max) - {lib_id}',
                cmap='viridis',
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_cadherin_minmax.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_cadherin_minmax.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 也可以绘制UMAP图查看整体分布
    sc.pl.umap(adata_spatial, 
               color=['hypoxia_adaptation_score_zscore'], 
               edges=False, 
               cmap='viridis')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.png", 
               dpi=300, bbox_inches='tight')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.pdf", 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print("缺氧适应基因得分计算完成并已生成图表")
    print(f"Z-score范围: {np.min(hypoxia_score_zscore):.3f} 到 {np.max(hypoxia_score_zscore):.3f}")
    print(f"Min-Max范围: {np.min(hypoxia_score_minmax):.3f} 到 {np.max(hypoxia_score_minmax):.3f}")
    
else:
    print("在数据中未找到任何缺氧适应相关基因")
# 为Z-score得分设置统一颜色尺度
hypoxia_zscore_min = adata_spatial.obs['hypoxia_adaptation_score_zscore'].min()
hypoxia_zscore_max = adata_spatial.obs['hypoxia_adaptation_score_zscore'].max()

# 为Min-Max得分设置统一颜色尺度（这个应该是0-1，但为了保险起见还是计算一下）
hypoxia_minmax_min = adata_spatial.obs['E2F1'].min()
hypoxia_minmax_max = adata_spatial.obs['hypoxia_adaptation_score_minmax'].max()
lib_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
# 绘制统一颜色尺度的缺氧适应基因得分图
for lib_id in lib_ids:
    try:
        ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
        
        # Z-score标准化结果
        sc.pl.spatial(
            ad,
            img_key="hires",
            color='E2F1',
            edges=False,
            title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
            cmap='viridis',
            library_id=lib_id,
            show=False,
            vmin=hypoxia_zscore_min,
            vmax=hypoxia_zscore_max
        )
        
        safe_lib_id = lib_id.replace('.', '_')
        plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_E2F1.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error processing hypoxia Z-score for {lib_id}: {str(e)}")
        plt.close()
        
cadherin_binding_up_genes = [
    'JUN'
]
available_genes = [gene for gene in cadherin_binding_up_genes if gene in adata_spatial.var_names]
print(f"在数据中找到的基因: {available_genes}")
print(f"缺失的基因: {set(cadherin_binding_up_genes) - set(available_genes)}")

# 提取这些基因的表达值
if len(available_genes) > 0:
    # 获取基因表达矩阵（使用原始counts或者已经标准化的数据）
    gene_expression = adata_spatial[:, available_genes].X
    
    # 转换为密集矩阵（如果需要）
    if not isinstance(gene_expression, np.ndarray):
        gene_expression = gene_expression.toarray()
    
    # 计算每个spot的基因表达总和
    hypoxia_score = np.sum(gene_expression, axis=1)
    
    # 标准化处理
    # 方法1: Z-score标准化（均值为0，标准差为1）
    hypoxia_score_zscore = (hypoxia_score - np.mean(hypoxia_score)) / np.std(hypoxia_score)
    
    # 方法2: Min-Max标准化（0-1范围）
    min_val = np.min(hypoxia_score)
    max_val = np.max(hypoxia_score)
    if max_val - min_val > 1e-10:
        hypoxia_score_minmax = (hypoxia_score - min_val) / (max_val - min_val)
    else:
        hypoxia_score_minmax = np.zeros_like(hypoxia_score)
    
    # 添加到obs中
    adata_spatial.obs['hypoxia_adaptation_score_zscore'] = hypoxia_score_zscore
    adata_spatial.obs['hypoxia_adaptation_score_minmax'] = hypoxia_score_minmax
    
    # 绘制每个样本的空间图
    library_ids = adata_spatial.obs['library_id'].unique()
    library_ids = ['ST_1.1.1', 'ST_1.2.1', 'ST_1.2.2', 'ST_1.3.1', 'ST_1.4.1']
    # 使用Z-score标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_zscore',
                edges=False,
                title=f'Hypoxia Adaptation Score (Z-score) - {lib_id}',
                cmap='viridis',  # 可以改为 'coolwarm', 'RdYlBu_r' 等其他colormap
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_E2F1.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_E2F1.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 使用Min-Max标准化的结果
    for lib_id in library_ids:
        try:
            # 过滤当前样本的数据
            ad = adata_spatial[adata_spatial.obs['library_id'] == lib_id].copy()
            
            # 绘制空间图
            sc.pl.spatial(
                ad,
                img_key="hires",
                color='hypoxia_adaptation_score_minmax',
                edges=False,
                title=f'Hypoxia Adaptation Score (Min-Max) - {lib_id}',
                cmap='viridis',
                library_id=lib_id,
                show=False
            )
            
            # 保存图像
            safe_lib_id = lib_id.replace('.', '_')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_E2F1_minmax.png", 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f"/home/data/sdzl14/NSCLC/zong/fig/spatial/{safe_lib_id}_E2F1_minmax.pdf", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error processing {lib_id}: {str(e)}")
            plt.close()
    
    # 也可以绘制UMAP图查看整体分布
    sc.pl.umap(adata_spatial, 
               color=['hypoxia_adaptation_score_zscore'], 
               edges=False, 
               cmap='viridis')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.png", 
               dpi=300, bbox_inches='tight')
    plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/spatial/umap_hypoxia_adaptation_zscore.pdf", 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print("缺氧适应基因得分计算完成并已生成图表")
    print(f"Z-score范围: {np.min(hypoxia_score_zscore):.3f} 到 {np.max(hypoxia_score_zscore):.3f}")
    print(f"Min-Max范围: {np.min(hypoxia_score_minmax):.3f} 到 {np.max(hypoxia_score_minmax):.3f}")
    
else:
    print("在数据中未找到任何缺氧适应相关基因")
