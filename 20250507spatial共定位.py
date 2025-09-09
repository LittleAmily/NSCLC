import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from turtle import title
from cytoolz import compose

import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss
# Set maximum number of jobs for Scanpy.
sc.settings.njobs = 32
RESOURCES_FOLDERNAME = "/home/data/sdzl14/NSCLC/zong/scenic/"
AUXILLIARIES_FOLDERNAME = "/home/data/sdzl14/NSCLC/zong/scenic/"
RESULTS_FOLDERNAME = "/home/data/sdzl14/NSCLC/zong/scenic/"
FIGURES_FOLDERNAME = "/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/scenic/"

BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
COLUMN_NAME_LOGO = "MotifLogo"
COLUMN_NAME_MOTIF_ID = "MotifID"
COLUMN_NAME_TARGETS = "TargetGenes"
def savesvg(fname: str, fig, folder: str=FIGURES_FOLDERNAME) -> None:

    """
    Save figure as vector-based SVG image format.
    """
    fig.tight_layout()
    fig.savefig(os.path.join(folder, fname), format='svg')
def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
    """
    :param df:
    :param base_url:
    """
    # Make sure the original dataframe is not altered.
    df = df.copy()
    
    # Add column with URLs to sequence logo.
    def create_url(motif_id):
        return '<img src=".png" style="max-height:124px;"></img>'.format(base_url, motif_id)
    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    
    # Truncate TargetGenes.
    def truncate(col_val):
        return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
    df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
    
    MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
    pd.set_option('display.max_colwidth', -1)
    display(HTML(df.head().to_html(escape=False)))
    pd.set_option('display.max_colwidth', MAX_COL_WIDTH)
import os, sys
os.getcwd()
os.listdir(os.getcwd())

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/malignant.h5ad')
row_attrs = {"Gene": np.array(x.var_names),}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create("sample.loom",x.X.transpose(),row_attrs,col_attrs)
# Downloaded fromm pySCENIC github repo: https://github.com/aertslab/pySCENIC/tree/master/resources
HUMAN_TFS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, 'lambert2018.txt')
# Ranking databases. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
RANKING_DBS_FNAMES = list(map(lambda fn: os.path.join(AUXILLIARIES_FOLDERNAME, fn),
                       ['hg19-500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather',
                       'hg19-tss-centered-5kb-10species.mc9nr.feather',
                        'hg19-tss-centered-10kb-10species.mc9nr.feather']))
# Motif annotations. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
MOTIF_ANNOTATIONS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, 'motifs-v9-nr.hgnc-m0.001-o0.0.tbl')
adata = x.copy()

METADATA_FNAME = os.path.join(RESULTS_FOLDERNAME, '.metadata.csv')
EXP_MTX_QC_FNAME = os.path.join(RESULTS_FOLDERNAME, '.qc.tpm.csv')
ADJACENCIES_FNAME = os.path.join(RESULTS_FOLDERNAME, '.adjacencies.tsv')
MOTIFS_FNAME = os.path.join(RESULTS_FOLDERNAME, '.motifs.csv')
REGULONS_DAT_FNAME = os.path.join(RESULTS_FOLDERNAME, '.regulons.dat')
AUCELL_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '.auc.csv')
BIN_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '.bin.csv')
THR_FNAME = os.path.join(RESULTS_FOLDERNAME, '.thresholds.csv')
ANNDATA_FNAME = os.path.join(RESULTS_FOLDERNAME, '.h5ad')
LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '.loom')
adata.write_h5ad(ANNDATA_FNAME) # Categorical dtypes are created.
adata.to_df().to_csv(EXP_MTX_QC_FNAME)
adata = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/adata_spatial.h5ad')
adata.obs['NK cells']
counts_file_path = "/home/data/sdzl14/NSCLC/zong/spatial/adata_spatial.h5ad"
adata = sc.read_h5ad(counts_file_path)
adata = adata.copy()
adata
p641 = adata[adata.obs['library_id'] == 'ST_6.4.1']
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import cdist

# 合并肿瘤细胞评分（假设评分为0-1的占比值）
p641.obs['tumor_total'] = p641.obs[['Cancer_P1/2/3', 'Cancer_P4', 'Cancer_P5']].sum(axis=1)

# 定义目标非肿瘤细胞类型列表
non_tumor_types = [
    'AT1', 'AT2', 'CD4+T', 'CD8+T', 'Ciliated cells', 
    'Club cells', 'Cycling cells', 'DC', 'Endothelial cells',
    'Epithelial cells', 'Fibroblasts', 'Macrophage', 'Mast cells',
    'Memory B cells', 'Monocytes', 'NK cells', 'Naive B cells',
    'Neutrophils', 'Plasma cells', 'pDC'
]

# 提取空间坐标（假设存储在obsm['spatial']）
coordinates = p641.obsm['spatial'][:, :2]  # 取前两列作为x,y坐标

# 计算所有spot之间的欧氏距离矩阵（优化内存版本）
def compute_distance_matrix(coords):
    """分块计算避免内存溢出"""
    chunk_size = 5000
    n = coords.shape[0]
    dist_mat = np.zeros((n, n))
    
    for i in range(0, n, chunk_size):
        for j in range(0, n, chunk_size):
            dist_mat[i:i+chunk_size, j:j+chunk_size] = cdist(
                coords[i:i+chunk_size], 
                coords[j:j+chunk_size], 
                'euclidean'
            )
    return dist_mat

distance_matrix = compute_distance_matrix(coordinates)
def calculate_weighted_distance(adata, cell_type, dist_matrix):
    """
    计算指定细胞类型与肿瘤细胞的加权平均距离
    
    参数:
        adata: AnnData对象
        cell_type: 目标非肿瘤细胞类型名称
        dist_matrix: 预计算的距离矩阵 (n_spot x n_spot)
    
    返回:
        (weighted_avg_distance, total_weight)
    """
    # 提取权重向量
    t = adata.obs['tumor_total'].values.reshape(-1, 1)  # (n_spot, 1)
    c = adata.obs[cell_type].values.reshape(1, -1)      # (1, n_spot)
    
    # 计算权重矩阵（外积）
    weights = t * c  # (n_spot, n_spot)
    
    # 计算加权距离总和（元素级乘积）
    weighted_sum = (weights * dist_matrix).sum()
    
    # 计算总权重
    total_weight = weights.sum()
    
    # 处理除零错误
    return weighted_sum / total_weight if total_weight > 0 else np.nan
# 存储结果的DataFrame
results = pd.DataFrame(index=non_tumor_types, columns=['weighted_distance', 'total_weight'])

# 逐细胞类型计算
for ct in non_tumor_types:
    try:
        dist = calculate_weighted_distance(p641, ct, distance_matrix)
        results.loc[ct, 'weighted_distance'] = dist
        print(f"Processed {ct}: {dist:.2f} μm")
    except KeyError:
        print(f"Warning: {ct} not found in adata.obs")
        continue

# 将结果存入AnnData的uns中
p641.uns['spatial_distance'] = results
# 结果可视化
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
results['weighted_distance'].sort_values().plot(kind='barh', color='steelblue')
plt.xlabel('Weighted Average Distance (μm)')
plt.title('Proximity Analysis: Tumor vs Other Cell Types')
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/ST_6.4.1.png', dpi=300)
##################  ###################   ############### ###############
# 获取所有library_id并排序
all_libs = adata.obs['library_id'].cat.categories.tolist()

# 初始化结果存储DataFrame (cell_types x library_ids)
all_results = pd.DataFrame(
    index=non_tumor_types,
    columns=all_libs,
    dtype=np.float32
)

# 主循环处理每个样本
for lib_id in all_libs:
    # 筛选当前样本
    p = adata[adata.obs['library_id'] == lib_id].copy()
    
    try:
        # 计算肿瘤总评分
        p.obs['tumor_total'] = p.obs[['Cancer_P1/2/3', 'Cancer_P4', 'Cancer_P5']].sum(axis=1)
        
        # 计算空间距离矩阵（优化内存）
        coords = p.obsm['spatial'][:, :2]
        distance_matrix = compute_distance_matrix(coords)
        
        # 并行计算所有细胞类型
        for ct in non_tumor_types:
            if ct not in p.obs.columns:  # 跳过不存在的列
                continue
                
            dist = calculate_weighted_distance(p, ct, distance_matrix)
            all_results.loc[ct, lib_id] = dist
            print(f"Processed {lib_id} - {ct}: {dist:.2f} μm")
            
        # 保存个体结果
        p.uns['spatial_distance'] = all_results[lib_id].copy()
        
    except Exception as e:
        print(f"Error processing {lib_id}: {str(e)}")
        continue

# 数据清洗：去除全NA的列/行
all_results.dropna(how='all', axis=0, inplace=True)
all_results.dropna(how='all', axis=1, inplace=True)

# 热图可视化
plt.figure(figsize=(12, 8))
sns.heatmap(
    all_results.T,  # 转置使样本为行，细胞类型为列
    cmap='coolwarm',
    annot=True,
    fmt=".1f",
    linewidths=.5,
    cbar_kws={'label': 'Weighted Distance (μm)'}
)
plt.title('Spatial Proximity: Tumor vs Cell Types Across Libraries')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/all_libs_heatmap.png', dpi=300, bbox_inches='tight')
# 数据归一化处理（0-10范围）
def scale_to_range(data, new_min=0, new_max=10):
    """ 最小-最大归一化 """
    data_min = data.min().min()  # 全局最小值
    data_max = data.max().max()  # 全局最大值
    
    # 处理全零或单一值情况
    if data_max == data_min:
        return pd.DataFrame(np.full_like(data, new_min), 
                          index=data.index, 
                          columns=data.columns)
    
    scaled = (data - data_min) / (data_max - data_min) * (new_max - new_min) + new_min
    return scaled

scaled_data = scale_to_range(all_results)

# 优化版热图绘制
plt.figure(figsize=(16, 12))
ax = sns.heatmap(
    scaled_data.T,
    cmap=sns.color_palette("rocket_r", as_cmap=True),  # 使用渐变红色系
    annot=True,
    annot_kws={
        'fontsize': 9,
        'fontweight': 'bold',
        'color': 'whitesmoke'  # 根据背景色调整字体颜色
    },
    fmt=".1f",
    linewidths=0.3,
    linecolor='lightgray',
    cbar_kws={
        'label': 'Normalized Distance (0-10)',
        'ticks': np.linspace(0, 10, 11)
    },
    vmin=0,
    vmax=10,
    square=True  # 使单元格为正方形
)

# 高级美化设置
ax.figure.axes[-1].yaxis.label.set_size(12)  # 颜色条标签字体
ax.figure.axes[-1].tick_params(labelsize=10)  # 颜色条刻度字体

# 坐标轴标签美化
ax.set_xticklabels(
    ax.get_xticklabels(),
    rotation=45,
    ha='right',
    fontstyle='italic',
    fontsize=11
)
ax.set_yticklabels(
    ax.get_yticklabels(), 
    fontsize=11,
    fontweight='semibold'
)

# 添加标题和说明
plt.title(
    "Normalized Spatial Proximity (0-10 Scale)\n",
    fontsize=14,
    fontweight='bold',
    pad=20
)
plt.xlabel(
    "Cell Type →",
    fontsize=12,
    labelpad=15,
    fontweight='bold'
)
plt.ylabel(
    "← Library ID",
    fontsize=12,
    labelpad=15,
    fontweight='bold'
)

# 添加网格参考线
ax.hlines(
    y=np.arange(scaled_data.shape[1])+0.5, 
    xmin=0, 
    xmax=scaled_data.shape[0], 
    colors='white',
    linewidth=0.1
)

# 保存结果
scaled_data.to_csv(
    '/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/normalized_0-10.csv'
)
plt.savefig(
    '/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/normalized_0-10_heatmap.png',
    dpi=350,
    bbox_inches='tight',
    facecolor='white'  # 确保白色背景
)
plt.savefig(
    '/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/normalized_0-10_heatmap.pdf',
    dpi=350,
    bbox_inches='tight',
    facecolor='white'  # 确保白色背景
)
# 解析library_id结构
def parse_library_id(lib_id):
    """ST_6.4.1 -> (patient=6, position=4, rep=1)"""
    parts = lib_id.split("_")[-1].split(".")
    return {
        'patient': int(parts[0]),
        'position': int(parts[1]),
        'replicate': int(parts[2])}

# 创建分层索引
meta_df = pd.DataFrame([parse_library_id(lib) for lib in all_libs], index=all_libs)

# 按病人统计（所有位置和重复的平均值）
patient_stats = scaled_data.T.groupby(meta_df['patient']).mean().T
patient_stats.columns = [f"Patient {p}" for p in patient_stats.columns]

# 按取材位置统计（所有病人和重复的平均值）
position_stats = scaled_data.T.groupby(meta_df['position']).mean().T
position_stats.columns = [f"Position {p}" for p in position_stats.columns]

# 可视化设置
def plot_stratified_heatmap(data, title, figsize=(14, 8)):
    plt.figure(figsize=figsize)
    ax = sns.heatmap(
        data,
        cmap='coolwarm',
        annot=True,
        fmt=".1f",
        linewidths=0.5,
        linecolor='lightgray',
        cbar_kws={'label': 'Normalized Distance'},
        annot_kws={'fontsize': 8}
    )
    
    # 坐标轴美化
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        ha='right',
        fontsize=10
    )
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0,
        fontsize=10
    )
    
    plt.title(title, fontsize=12, pad=20)
    plt.xlabel('')
    plt.ylabel('Cell Type', fontsize=10)
    return ax

# 绘制病人层面热图
plot_stratified_heatmap(
    patient_stats, 
    'Average Distance by Patient',
    figsize=(10, 12)
)
plt.savefig(
    '/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/patient_heatmap.png',
    dpi=300,
    bbox_inches='tight'
)
plt.savefig(
    '/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/patient_heatmap.pdf',
    dpi=300,
    bbox_inches='tight'
)

# 绘制取材位置层面热图
plot_stratified_heatmap(
    position_stats,
    'Average Distance by Sampling Position', 
    figsize=(8, 12)
)
plt.savefig(
    '/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/position_heatmap.png',
    dpi=300,
    bbox_inches='tight'
)
plt.savefig(
    '/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/position_heatmap.pdf',
    dpi=300,
    bbox_inches='tight'
)

# 可选：保存统计结果
patient_stats.to_csv('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/patient_stats.csv')
position_stats.to_csv('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/position_stats.csv')
# 创建癌症类型映射字典
cancer_type_map = {
    1: 'LUSC',
    2: 'LUSC',
    3: 'LUSC',
    5: 'LUAD',
    6: 'LUAD'
}

# 添加癌症类型元数据
meta_df['cancer_type'] = meta_df['patient'].map(cancer_type_map)

# 按癌症类型分组数据
lusc_data = scaled_data.loc[:, meta_df[meta_df['cancer_type'] == 'LUSC'].index]
luad_data = scaled_data.loc[:, meta_df[meta_df['cancer_type'] == 'LUAD'].index]

# 统计学分析函数
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

def calculate_significance(row):
    """计算两组间显著性"""
    lusc_vals = lusc_data.loc[row.name].dropna()
    luad_vals = luad_data.loc[row.name].dropna()
    
    if len(lusc_vals)<3 or len(luad_vals)<3:
        return np.nan
    
    # 使用Mann-Whitney U检验
    stat, p = mannwhitneyu(lusc_vals, luad_vals, alternative='two-sided')
    return p

# 计算所有细胞类型的p值
p_values = scaled_data.apply(calculate_significance, axis=1)
# FDR校正
reject, pvals_corrected, _, _ = multipletests(p_values, method='fdr_bh')
significance_df = pd.DataFrame({
    'p_raw': p_values,
    'p_adj': pvals_corrected,
    'significant': reject
})

# 计算各癌症类型均值
cancer_stats = pd.concat([
    lusc_data.mean(axis=1).rename('LUSC'),
    luad_data.mean(axis=1).rename('LUAD')
], axis=1).T

# 优化热图绘制
plt.figure(figsize=(16, 3))
ax = sns.heatmap(
    cancer_stats,
    cmap='viridis',
    annot=True,
    fmt=".1f",
    linewidths=0.5,
    linecolor='gray',
    cbar_kws={'label': 'Average Distance'},
    annot_kws={'fontsize':10, 'color':'white'},
    square=True
)

# 添加显著性标记
for i, ct in enumerate(cancer_stats.columns):
    if significance_df.loc[ct, 'significant']:
        ax.text(
            i + 0.5, 
            1.5,  # 在LUSC和LUAD行之间显示
            '*', 
            ha='center', 
            va='center', 
            color='red',
            fontsize=16,
            fontweight='bold'
        )

# 坐标轴美化
ax.set_xticklabels(
    ax.get_xticklabels(),
    rotation=45,
    ha='right',
    fontsize=10
)
ax.set_yticklabels(
    ax.get_yticklabels(),
    rotation=0,
    fontsize=12,
    fontweight='bold'
)

plt.title('Average Distance by Cancer Type\n* FDR < 0.05', pad=15)
plt.xlabel('')
plt.tight_layout()
plt.savefig(
    '/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/cancer_type_heatmap_dual.png',
    dpi=300,
    bbox_inches='tight'
)
plt.savefig(
    '/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/cancer_type_heatmap_dual.pdf',
    dpi=300,
    bbox_inches='tight'
)

# 绘制箱线图
melt_df = pd.concat([
    lusc_data.stack().reset_index().assign(cancer_type='LUSC'),
    luad_data.stack().reset_index().assign(cancer_type='LUAD')
], axis=0)
melt_df.columns = ['cell_type', 'sample', 'distance', 'cancer_type']

plt.figure(figsize=(15, 8))
ax = sns.boxplot(
    x='cell_type',
    y='distance',
    hue='cancer_type',
    data=melt_df,
    palette={'LUSC':'#1f77b4', 'LUAD':'#ff7f0e'},
    showfliers=False,
    width=0.7
)

# 添加统计标注
y_max = melt_df.groupby('cell_type')['distance'].max()
for i, ct in enumerate(scaled_data.index):
    if significance_df.loc[ct, 'significant']:
        ax.text(
            i, 
            y_max[ct] + 0.5, 
            f"p={significance_df.loc[ct, 'p_adj']:.2e}", 
            ha='center',
            color='red',
            fontsize=8
        )
        ax.plot([i-0.2, i+0.2], [y_max[ct]+0.3, y_max[ct]+0.3], color='red', lw=1)

plt.title('Distance Distribution by Cancer Type', pad=20)
plt.xlabel('Cell Type')
plt.ylabel('Normalized Distance')
plt.xticks(rotation=45)
plt.legend(title='Cancer Type', loc='upper right')
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/cancer_type_boxplot.png', dpi=300)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/cancer_type_boxplot.pdf', dpi=300)
# 保存统计结果
significance_df.to_csv('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/cancer_type_stats.csv')
# 按细胞类型相似性排序
g = sns.clustermap(
    cancer_stats,
    row_cluster=False,
    col_cluster=True,
    figsize=(16, 3),
    cmap='viridis'
)
g.ax_heatmap.set_title('Clustered by Cell Type Similarity')
g.savefig(
    '/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/cancer_type_heatmap_clustered.png',
    dpi=300,
    bbox_inches='tight'
)
g.savefig(
    '/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/cancer_type_heatmap_clustered.pdf',
    dpi=300,
    bbox_inches='tight'
)
# 定义position映射
position_map = {
    1: 'tumor_edge',
    2: 'tumor_middle',
    3: 'normal_adjacent',
    4: 'normal_distance'
}
meta_df['position_name'] = meta_df['position'].map(position_map)

# 按position分组数据
position_groups = all_results.T.groupby(meta_df['position_name'])

# 统计检验
from scipy.stats import kruskal
from scikit_posthocs import posthoc_dunn

position_stats = pd.DataFrame(index=all_results.index, columns=['p_value', 'significance'])
# 可视化
def format_significance(pval_matrix):
    """将Dunn检验结果格式化为符号表示"""
    symbols = []
    for (i, j) in combinations(pval_matrix.index, 2):
        if pval_matrix.loc[i, j] < 0.001:
            symbols.append(f"{i} vs {j} ***")
        elif pval_matrix.loc[i, j] < 0.01:
            symbols.append(f"{i} vs {j} **")
        elif pval_matrix.loc[i, j] < 0.05:
            symbols.append(f"{i} vs {j} *")
    return '\n'.join(symbols)
for ct in scaled_data.index:
    group_data = [group[ct].dropna() for name, group in position_groups if ct in group]
    
    # 至少需要2组且有数据才进行检验
    if len(group_data) < 2 or sum(len(g)>=3 for g in group_data) < 2:
        position_stats.loc[ct] = [np.nan, '']
        continue
    
    # Kruskal-Wallis H检验
    h_stat, p_val = kruskal(*group_data)
    position_stats.loc[ct, 'p_value'] = p_val
    
    # 事后Dunn检验（多重校正）
    if p_val < 0.05:
        melted = scaled_data.loc[[ct]].T.join(meta_df['position_name']).dropna()
        dunn_pvals = posthoc_dunn(melted, val_col=ct, group_col='position_name', p_adjust='fdr_bh')
        position_stats.loc[ct, 'significance'] = format_significance(dunn_pvals)


# 热图绘制
position_means = scaled_data.T.groupby(meta_df['position_name']).mean().T
order = ['tumor_edge', 'tumor_middle', 'normal_adjacent', 'normal_distance']

plt.figure(figsize=(12, 8))
ax = sns.heatmap(
    position_means[order],
    cmap='YlGnBu',
    annot=True,
    fmt=".1f",
    linewidths=0.5,
    cbar_kws={'label': 'Average Distance'},
    annot_kws={'fontsize':9}
)

# 标注显著性
for idx, ct in enumerate(position_means.index):
    if position_stats.loc[ct, 'p_value'] < 0.05:
        ax.text(4.2, idx+0.5, position_stats.loc[ct, 'significance'], 
                va='center', fontsize=8, color='maroon')

plt.title('Spatial Distance Across Sampling Positions\n(右侧标注组间差异显著性)', pad=20)
plt.xlabel('Sampling Position')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/position_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/position_heatmap.pdf', dpi=300, bbox_inches='tight')
# 箱线图绘制
melt_df = scaled_data.stack().reset_index()
melt_df.columns = ['cell_type', 'sample', 'distance']
melt_df = melt_df.join(meta_df['position_name'], on='sample')

plt.figure(figsize=(18, 10))
sns.boxplot(
    x='cell_type',
    y='distance',
    hue='position_name',
    data=melt_df,
    hue_order=order,
    palette='Set2',
    showfliers=False
)

# 自动标注显著性
y_max = melt_df.groupby('cell_type')['distance'].max()
for idx, ct in enumerate(position_means.index):
    if position_stats.loc[ct, 'p_value'] < 0.05:
        plt.text(
            idx, 
            y_max[ct]+0.8, 
            '*', 
            ha='center', 
            color='red',
            fontsize=14
        )

plt.title('Distance Distribution by Sampling Position', pad=20)
plt.xlabel('Cell Type')
plt.ylabel('Normalized Distance')
plt.xticks(rotation=45)
plt.legend(title='Position', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/position_boxplot.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/position_boxplot.pdf', dpi=300, bbox_inches='tight')

# 保存结果
position_stats.to_csv('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/position_statistics.csv')
# 筛选tumor_edge样本
tumor_edge_samples = meta_df[meta_df['position_name'] == 'tumor_middle'].index
te_data = all_results[tumor_edge_samples]

# 创建癌症类型映射
te_meta = meta_df.loc[tumor_edge_samples, 'cancer_type']

# 数据归一化 (1-10范围)
def scale_1_10(data):
    data_min = data.min().min()
    data_max = data.max().max()
    if data_max == data_min:
        return pd.DataFrame(np.full_like(data, 1), index=data.index, columns=data.columns)
    return (data - data_min)/(data_max - data_min)*9 + 1  # 1-10 scaling

scaled_te = scale_1_10(te_data)

# 统计学检验
from scipy.stats import mannwhitneyu
te_stats = pd.DataFrame(index=scaled_te.index, columns=['LUSC_median', 'LUAD_median', 'p_value', 'adj_p'])

for ct in scaled_te.index:
    lusc = scaled_te.loc[ct, te_meta == 'LUSC'].dropna()
    luad = scaled_te.loc[ct, te_meta == 'LUAD'].dropna()
    
    if len(lusc)<3 or len(luad)<3:
        te_stats.loc[ct] = [np.nan]*4
        continue
    
    te_stats.loc[ct, 'LUSC_median'] = lusc.median()
    te_stats.loc[ct, 'LUAD_median'] = luad.median()
    _, p = mannwhitneyu(lusc, luad)
    te_stats.loc[ct, 'p_value'] = p

# FDR校正
te_stats['adj_p'] = multipletests(te_stats['p_value'], method='fdr_bh')[1]

# 热图绘制
plot_data = pd.concat([
    scaled_te.loc[:, te_meta == 'LUSC'].mean(axis=1).rename('LUSC'),
    scaled_te.loc[:, te_meta == 'LUAD'].mean(axis=1).rename('LUAD')
], axis=1)

plt.figure(figsize=(8, 10))
ax = sns.heatmap(
    plot_data,
    cmap='coolwarm',
    annot=True,
    fmt=".1f",
    linewidths=0.5,
    cbar_kws={'label': 'Scaled Distance (1-10)'},
    annot_kws={'fontsize':9, 'fontweight':'bold'},
    mask=plot_data.isnull()
)

# 标注显著性
for idx, ct in enumerate(plot_data.index):
    if te_stats.loc[ct, 'adj_p'] < 0.05:
        ax.text(2.5, idx+0.5, f"*", 
                ha='center', va='center', 
                color='black', fontsize=14, fontweight='bold')

plt.title('Tumor Edge: LUSC vs LUAD Comparison\n* FDR < 0.05', pad=15)
plt.xlabel('Cancer Type')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/tumor_edge_heatmap.png', 
           dpi=300, bbox_inches='tight')

# 箱线图绘制
melt_te = scaled_te.stack().reset_index()
melt_te.columns = ['cell_type', 'sample', 'distance']
melt_te = melt_te.join(te_meta, on='sample')

plt.figure(figsize=(15, 8))
ax = sns.boxplot(
    x='cell_type',
    y='distance',
    hue='cancer_type',
    data=melt_te,
    palette={'LUSC':'#1f77b4', 'LUAD':'#ff7f0e'},
    showfliers=False,
    width=0.7,
    linewidth=1.5
)

# 显著性标注
y_max = melt_te.groupby('cell_type')['distance'].max()
for idx, ct in enumerate(plot_data.index):
    if te_stats.loc[ct, 'adj_p'] < 0.05:
        ax.text(
            idx, 
            y_max[ct] + 0.5, 
            f"p={te_stats.loc[ct, 'adj_p']:.2e}", 
            ha='center',
            color='red',
            fontsize=8
        )
        ax.plot([idx-0.2, idx+0.2], [y_max[ct]+0.3, y_max[ct]+0.3], 
                color='red', lw=1.5)

plt.title('Tumor Edge Spatial Distance Distribution', pad=20)
plt.xlabel('Cell Type')
plt.ylabel('Scaled Distance (1-10)')
plt.xticks(rotation=45)
plt.legend(title='Cancer Type', loc='upper right')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/tumor_edge_boxplot.png', 
          dpi=300, bbox_inches='tight')

# 保存结果
te_stats.to_csv('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/tumor_edge_stats.csv')
# 计算差异方向
te_stats['delta'] = te_stats['LUAD_median'] - te_stats['LUSC_median']
te_stats['trend'] = te_stats['delta'].apply(lambda x: 'LUAD > LUSC' if x>0 else 'LUSC > LUAD')

# 绘制差异方向分布
plt.figure(figsize=(8,4))
sns.countplot(data=te_stats, x='trend', order=['LUAD > LUSC', 'LUSC > LUAD'])
plt.title('Differential Trend in Tumor Edge')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/tumor_edge_significant_boxplot.png', 
          dpi=300, bbox_inches='tight')
# 筛选LUSC患者数据
lusc_samples = meta_df[meta_df['cancer_type'] == 'LUSC'].index
lusc_data = all_results[lusc_samples]

# 独立归一化 (1-10范围)
def scale_1_10_lusc(data):
    data_min = data.min().min()
    data_max = data.max().max()
    if data_max == data_min:
        return pd.DataFrame(np.full_like(data, 1), index=data.index, columns=data.columns)
    return (data - data_min)/(data_max - data_min)*9 + 1

scaled_lusc = scale_1_10_lusc(lusc_data)

# 统计检验
from scipy.stats import kruskal
from scikit_posthocs import posthoc_dunn

position_order = ['tumor_edge', 'tumor_middle', 'normal_adjacent', 'normal_distance']
lusc_stats = pd.DataFrame(index=scaled_lusc.index, columns=['p_value', 'significance'])

for ct in scaled_lusc.index:
    group_data = []
    for pos in position_order:
        samples = meta_df[(meta_df['position_name'] == pos) & (meta_df.index.isin(lusc_samples))].index
        group = scaled_lusc.loc[ct, samples].dropna()
        if len(group) >= 3:  # 至少3个样本才纳入分析
            group_data.append(group)
    
    if len(group_data) < 2:
        lusc_stats.loc[ct] = [np.nan, '']
        continue
    
    # Kruskal-Wallis检验
    h_stat, p_val = kruskal(*group_data)
    lusc_stats.loc[ct, 'p_value'] = p_val
    
    # 事后检验
    if p_val < 0.05:
        melted = scaled_lusc.loc[[ct]].T.join(meta_df['position_name']).dropna()
        dunn_pvals = posthoc_dunn(melted, val_col=ct, group_col='position_name', p_adjust='fdr_bh')
        # 格式化显著性标记
        sig_symbols = []
        for (i, j) in combinations(dunn_pvals.index, 2):
            p = dunn_pvals.loc[i, j]
            if p < 0.001:
                sig_symbols.append(f"{i[:2]}vs{j[:2]}***")  # 简写位置名称
            elif p < 0.01:
                sig_symbols.append(f"{i[:2]}vs{j[:2]}**")
            elif p < 0.05:
                sig_symbols.append(f"{i[:2]}vs{j[:2]}*")
        lusc_stats.loc[ct, 'significance'] = ', '.join(sig_symbols)

# 热图绘制
position_means = scaled_lusc.T.groupby(meta_df['position_name']).mean().T[position_order]

plt.figure(figsize=(12, 8))
ax = sns.heatmap(
    position_means,
    cmap='YlOrRd',
    annot=True,
    fmt=".1f",
    linewidths=0.5,
    cbar_kws={'label': 'Scaled Distance (1-10)'},
    annot_kws={'fontsize':9, 'fontweight':'bold'},
    mask=position_means.isnull()
)

# 标注显著性
for idx, ct in enumerate(position_means.index):
    if pd.notna(lusc_stats.loc[ct, 'p_value']) and lusc_stats.loc[ct, 'p_value'] < 0.05:
        ax.text(
            4.5, idx+0.5, 
            lusc_stats.loc[ct, 'significance'],
            va='center',
            fontsize=8,
            color='darkblue',
            fontstyle='italic'
        )

plt.title('LUSC: Spatial Distance Across Sampling Positions\n(Right: Significant Comparisons)', pad=20)
plt.xlabel('Sampling Position')
plt.ylabel('Cell Type')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/LUSC_position_heatmap.png', 
           dpi=300, bbox_inches='tight')

# 箱线图绘制
melt_lusc = scaled_lusc.stack().reset_index()
melt_lusc.columns = ['cell_type', 'sample', 'distance']
melt_lusc = melt_lusc.join(meta_df['position_name'], on='sample')

plt.figure(figsize=(18, 10))
ax = sns.boxplot(
    x='cell_type',
    y='distance',
    hue='position_name',
    data=melt_lusc,
    hue_order=position_order,
    palette='Paired',
    showfliers=False,
    linewidth=1.2
)

# 显著性标注
y_max = melt_lusc.groupby('cell_type')['distance'].max()
for idx, ct in enumerate(position_means.index):
    if pd.notna(lusc_stats.loc[ct, 'p_value']) and lusc_stats.loc[ct, 'p_value'] < 0.05:
        ax.text(
            idx, 
            y_max[ct] + 0.8, 
            '*', 
            ha='center', 
            color='red',
            fontsize=14,
            fontweight='bold'
        )

plt.title('LUSC: Distance Distribution by Sampling Position', pad=20)
plt.xlabel('Cell Type')
plt.ylabel('Scaled Distance (1-10)')
plt.xticks(rotation=45)
plt.legend(title='Position', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/LUSC_position_boxplot.png', 
          dpi=300, bbox_inches='tight')

# 保存结果
lusc_stats.to_csv('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/LUSC_position_stats.csv')
def calculate_normalized_distance(adata, cell_type, dist_matrix):
    """
    计算标准化后的肿瘤-非肿瘤细胞加权距离
    包含三项改进：
    1. 肿瘤评分合并
    2. 权重归一化
    3. 距离密度校正
    """
    # 合并肿瘤细胞评分 (所有肿瘤亚型合并)
    adata.obs['tumor_combined'] = adata.obs[['Cancer_P1/2/3', 'Cancer_P4', 'Cancer_P5']].sum(axis=1)
    
    # 提取权重向量并进行归一化
    tumor_scores = adata.obs['tumor_combined'].values
    cell_scores = adata.obs[cell_type].values
    
    # 密度校正因子（避免spot数量影响）
    density_factor = np.sqrt(np.mean(tumor_scores)) * np.sqrt(np.mean(cell_scores))
    
    # 计算权重矩阵（外积）并进行双向归一化
    weights = np.outer(tumor_scores, cell_scores)
    weights /= (weights.sum(axis=1, keepdims=True) + 1e-8)  # 行归一化
    weights /= (weights.sum(axis=0, keepdims=True) + 1e-8)  # 列归一化
    
    # 计算校正后的加权距离
    weighted_dist = (weights * dist_matrix).sum()
    normalized_dist = weighted_dist / (density_factor + 1e-8)
    
    return normalized_dist
# 预计算距离矩阵（优化内存）
coordinates = p641.obsm['spatial'][:, :2]
distance_matrix = compute_distance_matrix(coordinates)

# 定义非肿瘤细胞类型列表
non_tumor_types = [
    'AT1', 'AT2', 'CD4+T', 'CD8+T', 'Ciliated cells',
    'Club cells', 'Cycling cells', 'DC', 'Endothelial cells',
    'Epithelial cells', 'Fibroblasts', 'Macrophage', 'Mast cells',
    'Memory B cells', 'Monocytes', 'NK cells', 'Naive B cells',
    'Neutrophils', 'Plasma cells', 'pDC'
]

# 并行计算（使用Joblib加速）
from joblib import Parallel, delayed

results = pd.DataFrame(index=non_tumor_types, columns=['normalized_distance'])

def process_cell_type(ct):
    try:
        if ct not in p641.obs.columns:
            print(f"Warning: {ct} not found, skipping...")
            return (ct, np.nan)
            
        dist = calculate_normalized_distance(p641, ct, distance_matrix)
        print(f"Processed {ct}: {dist:.2f} (Normalized Units)")
        return (ct, dist)
    except Exception as e:
        print(f"Error processing {ct}: {str(e)}")
        return (ct, np.nan)

# 使用4个CPU核心并行计算
outputs = Parallel(n_jobs=4)(delayed(process_cell_type)(ct) for ct in non_tumor_types)

# 整理结果
for ct, dist in outputs:
    results.loc[ct, 'normalized_distance'] = dist

# 保存结果
p641.uns['spatial_distance_normalized'] = results
# 绘制距离与细胞密度相关性
densities = p641.obs[non_tumor_types].mean()
distances = results['weighted_distance']

# 确保 distances 是数值类型
distances = pd.to_numeric(distances, errors='coerce')

# 创建一个新的 DataFrame 包含 densities 和 distances
df = pd.DataFrame({'density': densities, 'distance': distances}).dropna()

# 进行线性回归分析
X = df['density']
y = df['distance']

# 添加常数项
import statsmodels.api as sm
X = sm.add_constant(X)

# 拟合模型
model = sm.OLS(y, X).fit()

# 获取 R² 和 p 值
r_squared = model.rsquared
p_value = model.pvalues[1]

# 绘制回归图
plt.figure(figsize=(10, 6))
sns.regplot(x=df['density'], y=df['distance'], ci=95, scatter_kws={'alpha':0.7})
plt.xlabel('Density')
plt.ylabel('Distance')
plt.title(f'Regression Plot of Density vs Distance\n$R^2$: {r_squared:.2f}, P-value: {p_value:.2g}')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/ST_6.4.1_regplot.png', dpi=300)

# 结果可视化
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
results['normalized_distance'].sort_values().plot(kind='barh', color='steelblue')
plt.xlabel('Weighted Average Distance (μm)')
plt.title('Proximity Analysis: Tumor vs Other Cell Types')
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/ST_6.4.1_regress.png', dpi=300)
import squidpy as sq
sq.pl.spatial(p641, color=['tumor_combined', 'Cancer_P1/2/3', 'Cancer_P4', 'Cancer_P5'],
              size=0.5, alpha_img=0.5, alpha_overlay=0.5,
              save='/home/data/sdzl14/NSCLC/zong/fig/spatial/ST_6.4.1_spatial.png')
import matplotlib.pyplot as plt



# 设置可视化参数
sc.set_figure_params(figsize=(12, 6), dpi=300)
fig, axs = plt.subplots(ncols=2, gridspec_kw={'wspace':0.3})
p311 = adata[adata.obs['library_id'] == 'ST_3.1.1']
# 计算距离梯度指标
# 合并肿瘤细胞评分并创建可视化字段
p311.obs['tumor_combined'] = p311.obs[['Cancer_P1/2/3', 'Cancer_P4', 'Cancer_P5']].sum(axis=1)
p311.obs['Macrophage_tumor'] = p311.obs['Macrophage'] * p311.obs['tumor_combined']
p311.obs['CD8T_tumor'] = p311.obs['CD8+T'] * p311.obs['tumor_combined']

# 设置可视化参数
sc.set_figure_params(figsize=(12, 6), dpi=300)
fig, axs = plt.subplots(ncols=2, gridspec_kw={'wspace':0.3})

# 绘制Macrophage与肿瘤的共定位
sc.pl.spatial(
    p311,
    color='Macrophage_tumor',
    alpha_img=0.3,
    size=1.5,
    color_map='Blues',
    vmax='p99',  # 使用99分位数作为最大值
    title='Macrophage-Tumor Proximity',
    ax=axs[0],
    show=False,
    library_id='ST_3.1.1'
)
# 叠加肿瘤细胞分布
sc.pl.spatial(
    p311,
    color='tumor_combined',
    alpha=0.5,  # 降低透明度避免遮挡
    size=0.8,
    color_map='Greens',
    ax=axs[0],
    show=False,
    library_id='ST_3.1.1',
    title='Macrophage-Tumor Proximity',
)

# 绘制CD8+T与肿瘤的共定位
sc.pl.spatial(
    p311,
    color='CD8T_tumor',
    alpha_img=0.3,
    size=1.5,
    color_map='Blues',
    vmax='p99',
    title='CD8+T-Tumor Proximity',
    ax=axs[1],
    show=False,
    library_id='ST_3.1.1'
)
# 叠加肿瘤细胞分布
sc.pl.spatial(
    p311,
    color='tumor_combined',
    alpha=0.5,
    size=0.8,
    color_map='Greens',
    ax=axs[1],
    show=False,
    library_id='ST_3.1.1',
    title='CD8+T-Tumor Proximity',
)

# 添加自定义图例
for ax in axs:
    ax.scatter([], [], c='blue', s=30, label='High Co-localization')
    ax.scatter([], [], c='green', s=15, label='Tumor Cells')
    ax.legend(loc='upper left', frameon=False)

plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/spatial/ST_3.1.1_spatial_comparison.png', bbox_inches='tight')
import matplotlib.pyplot as plt
from math import ceil

# 配置参数
non_tumor_types = [
    'CD4+T', 'CD8+T', 
     'DC', 'Endothelial cells',
    'Epithelial cells', 'Fibroblasts', 'Macrophage', 'Mast cells',
    'Memory B cells', 'Monocytes', 'NK cells', 'Naive B cells',
    'Neutrophils', 'Plasma cells','pDC']
all_libs = adata.obs['library_id'].unique().tolist()
n_cols = 4  # 每行显示4个细胞类型
color_params = {
    'co_localization': {'cmap': 'Blues', 'vmax': 'p99'},
    'tumor': {'cmap': 'Greys', 'alpha': 0.3, 'size': 0.8}
}
# 预计算全局颜色范围
def calculate_global_vmax(adata, non_tumor_types):
    """计算所有样本所有细胞类型的99分位数最大值"""
    global_max = 0
    for lib_id in adata.obs['library_id'].unique():
        p = adata[adata.obs['library_id'] == lib_id].copy()
        p.obs['tumor_combined'] = p.obs[['Cancer_P1/2/3', 'Cancer_P4', 'Cancer_P5']].sum(axis=1)
        
        for ct in non_tumor_types:
            if ct not in p.obs.columns:
                continue
            col_name = f"{ct}_tumor".replace('+', '_')
            p.obs[col_name] = p.obs[ct] * p.obs['tumor_combined']
            
            # 计算当前细胞类型的99分位数
            current_max = p.obs[col_name].quantile(0.99)
            global_max = max(global_max, current_max)
    
    return global_max * 1.05  # 增加5%缓冲空间

# 计算全局vmax
global_vmax = calculate_global_vmax(adata, non_tumor_types)
global_vmax = 100
# 更新颜色参数
color_params = {
    'co_localization': {
        'cmap': 'viridis',  # 改用更通用的色阶
        'vmin': 0,
        'vmax': global_vmax,
        'colorbar_loc': None  # 后续统一添加颜色条
    },
    'tumor': {'cmap': 'Greys', 'alpha': 0.3, 'size': 0.8}
}
# 修正后的主循环部分
for lib_id in all_libs:
    p = adata[adata.obs['library_id'] == lib_id].copy()
    if p.n_obs == 0:
        print(f"跳过不存在的样本: {lib_id}")
        continue
    
    # 预处理
    p.obs['tumor_combined'] = p.obs[['Cancer_P1/2/3', 'Cancer_P4', 'Cancer_P5']].sum(axis=1)
    
    # 动态计算子图布局
    n_types = len(non_tumor_types)
    n_rows = ceil(n_types / n_cols)
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols*5, n_rows*4))
    axs = axs.ravel()
    
    # 添加颜色条
    cax = fig.add_axes([0.92, 0.3, 0.02, 0.4])
    norm = mpl.colors.Normalize(vmin=0, vmax=global_vmax)
    cmap = plt.get_cmap(color_params['co_localization']['cmap'])
    mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical', label='Co-localization Score')
    
    for i, ct in enumerate(non_tumor_types):
        ax = axs[i]
        col_name = f"{ct}_tumor".replace('+', '_')  # 先定义col_name
        
        try:
            # 创建共定位评分
            p.obs[col_name] = p.obs[ct] * p.obs['tumor_combined']
            
            # 绘制共定位图
            sc.pl.spatial(
                p,
                color=col_name,
                ax=ax,
                show=False,
                title=f"{ct}\n(max: {p.obs[col_name].max():.1f})",  # 在标题中显示最大值
                size=1.5,
                alpha_img=0.3,
                library_id=lib_id,
                **color_params['co_localization']
            )
            
            # 叠加肿瘤分布
            sc.pl.spatial(
                p,
                color='tumor_combined',
                ax=ax,
                show=False,
                **color_params['tumor'],
                library_id=lib_id,
                title=f"{ct}\n(max: {p.obs['tumor_combined'].max():.1f})",
            )
            
            # 添加图例
            ax.scatter([], [], c='royalblue', s=30, label='Co-localization')
            ax.scatter([], [], c='dimgrey', s=15, label='Tumor')
            ax.legend(loc='upper right', frameon=False, fontsize=8)
            
        except KeyError as e:
            print(f"警告：{ct} 或相关字段在样本 {lib_id} 中不存在，跳过绘制")
            ax.axis('off')
            continue
    
    # 隐藏空子图
    for j in range(i+1, len(axs)):
        axs[j].axis('off')
    
    plt.suptitle(f"Sample {lib_id} Tumor Co-localization Patterns", y=1.02, fontsize=14)
    plt.tight_layout()
    
    # 保存
    save_path = f'/home/data/sdzl14/NSCLC/zong/fig/spatial/distance/{lib_id}_co_localization.pdf'
    plt.savefig(save_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"已保存: {save_path}")