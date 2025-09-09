from email import header
from turtle import color
from click import group
from matplotlib import legend
from pyparsing import col
import scanpy as sc

import matplotlib.pyplot as plt
import scvelo as scv
import cellrank as cr
stromal = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/stromal.h5ad')
stromal = stromal.copy()
stromal.obs['stromal_celltype'].cat.categories
p = ['ADH1B+ CAF', 'Artery', 'BCHE+ SMC', 'COL11A1+ CAF',
       'CPE+ Venule', 'Capillary', 'Lymphatic EC', 'MYH11+ Pericyte',
       'Pericyte', 'SELE+ Venule', 'SMC', 'Tip', 'Venule']
stromal =  stromal[stromal.obs['stromal_celltype'].isin(p)]
import palantir
sc.pp.normalize_total(stromal, target_sum=1e4)
sc.pp.log1p(stromal)
sc.pp.highly_variable_genes(stromal, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(stromal, max_value=10)
sc.pp.neighbors(stromal, use_rep='X_scanvi_fix_linear',n_neighbors=35)
sc.tl.umap(stromal,min_dist=0.5)
sc.pl.umap(stromal,color=['stromal_celltype'])
sc.tl.diffmap(stromal)
stromal.obsm['X_diffmap_'] = stromal.obsm['X_diffmap'][:,1:]
sc.pl.embedding(stromal,'diffmap_',color=['stromal_celltype'])
results_file = '/home/data/sdzl14/NSCLC/zong/fig/stromal/'
plt.savefig(results_file + 'diffmap_stromal_celltype-.png', dpi=300, bbox_inches='tight')
stromal.obsm['X_diffmap'].shape
#(65680, 15)
sc.tl.pca(stromal, svd_solver='arpack')
sc.pp.neighbors(stromal, n_neighbors=15, n_pcs=15, use_rep='X_diffmap')
sc.tl.draw_graph(stromal)
sc.pl.draw_graph(stromal, color=['stromal_celltype'])
sc.tl.diffmap(stromal)
sc.pp.neighbors(stromal, n_neighbors=20, use_rep='X_diffmap_')
sc.tl.draw_graph(stromal)
sc.pl.draw_graph(stromal, color='stromal_celltype')

plt.savefig(results_file + 'diffmap——draw_graph_stromal_celltype.png', dpi=300, bbox_inches='tight')

# 可视化
sc.pl.embedding(stromal, basis="force_directed")

sc.pp.neighbors(stromal, n_neighbors=10, n_pcs=10, use_rep='X_diffmap_')
sc.tl.umap(stromal)
sc.pl.umap(stromal,color='stromal_celltype')
plt.savefig(results_file + 'umap_stromal_celltype.png', dpi=300, bbox_inches='tight')
sc.tl.leiden(stromal, resolution=0.3)
sc.tl.paga(stromal, groups='leiden')
sc.pl.paga(stromal, color='leiden')
plt.savefig(results_file + 'paga_leiden.png', dpi=300, bbox_inches='tight')
sc.tl.paga(stromal, groups='stromal_celltype')
sc.pl.paga(stromal, color='stromal_celltype')
plt.savefig(results_file + 'paga_stromal_celltype.png', dpi=300, bbox_inches='tight')
sc.pl.umap(stromal,color='leiden')
plt.savefig(results_file + 'umap_leiden.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'umap_leiden.pdf', dpi=300, bbox_inches='tight')
stromal.write_h5ad('/home/data/sdzl14/NSCLC/zong/stromal.h5ad')
import pandas as pd
df = pd.read_csv('/home/data/sdzl14/NSCLC/zong/stromal.csv')
print(df.columns)
# 提取需要保存的列
columns_to_save = ['CytoTRACE2_Score', 'CytoTRACE2_Potency', 'CytoTRACE2_Relative', 
                   'preKNN_CytoTRACE2_Score', 'preKNN_CytoTRACE2_Potency']

# 将这些列添加到 stromal.obs 中
stromal.obs[columns_to_save] = df.set_index('Unnamed: 0').loc[stromal.obs.index, columns_to_save]

sc.pl.umap(stromal, color='CytoTRACE2_Score', cmap='plasma')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_CytoTRACE2_Score.png', dpi=300, bbox_inches='tight')
sc.pl.umap(stromal, color='preKNN_CytoTRACE2_Score', cmap='plasma')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_preKNN_CytoTRACE2_Score.png', dpi=300, bbox_inches='tight')
sc.pl.umap(stromal, color='stromal_celltype')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_stromal_celltype.png', dpi=300, bbox_inches='tight')

sc.pl.umap(stromal, color='preKNN_CytoTRACE2_Potency')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_preKNN_CytoTRACE2_Potency.png', dpi=300, bbox_inches='tight')
sc.pl.umap(stromal, color='GPNMB')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_GPNMB.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_GPNMB.pdf', dpi=300, bbox_inches='tight')
sc.pl.umap(stromal, color='SET')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_SET.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/umap_stromal_celltype.png', dpi=300, bbox_inches='tight')
# 提取所需列
df = stromal.obs[['stromal_celltype', 'preKNN_CytoTRACE2_Score']]

# 计算每个 celltype 组内的 preKNN_CytoTRACE2_Score 均值，并按均值排序
order = df.groupby('stromal_celltype')['preKNN_CytoTRACE2_Score'].median().sort_values(ascending=False).index
import seaborn as sns
# 设置配色
palette = sns.color_palette("Spectral", len(order))

plt.figure(figsize=(12, 6))

# 小提琴图（按 order 排序）
sns.violinplot(x='stromal_celltype', y='preKNN_CytoTRACE2_Score', data=df,
               order=order, palette=palette, inner=None, linewidth=1.2)

# 叠加箱型图（按 order 排序）
sns.boxplot(x='stromal_celltype', y='preKNN_CytoTRACE2_Score', data=df,
            order=order,
            showcaps=False,
            boxprops={'facecolor': 'none'},
            showfliers=False,
            whiskerprops={'linewidth': 1.5},
            width=0.4)

# 图像美化
plt.xticks(rotation=45, fontsize=10)
plt.title('preKNN_CytoTRACE2_Score by Stromal Celltype (Ordered by Median)', fontsize=14)
plt.ylabel('preKNN_CytoTRACE2_Score', fontsize=12)
plt.xlabel('Stromal Celltype', fontsize=12)
plt.tight_layout()

# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/stromal/'
plt.savefig(results_file + 'violin_boxplot_preKNN_CytoTRACE2_Score_sorted.png',
            dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'violin_boxplot_preKNN_CytoTRACE2_Score_sorted.pdf',
            dpi=300, bbox_inches='tight')
import numpy as np
stromal.uns["iroot"] = np.flatnonzero(stromal.obs["stromal_celltype"] == "SMC")[0]

sc.tl.dpt(stromal)
sc.pl.draw_graph(stromal, color=["dpt_pseudotime"], legend_loc="on data")
plt.savefig(results_file + 'stromal_dpt_pseudotime.png',dpi = 300, bbox_inches='tight')
plt.savefig(results_file + 'stromal_dpt_pseudotime.pdf',dpi = 300, bbox_inches='tight')
sc.pl.paga(stromal,color=["dpt_pseudotime"])
plt.savefig(results_file + 'stromal_dpt_pseudotime_paga.png',dpi = 300, bbox_inches='tight')
plt.savefig(results_file + 'stromal_dpt_pseudotime_paga.pdf',dpi = 300, bbox_inches='tight')
stromal.obs['stromal_celltype'].cat.categories
caf = ['ADH1B+ CAF','COL11A1+ CAF']
CAF = stromal[stromal.obs['stromal_celltype'].isin(caf)]
CAF.X = CAF.layers['counts']

sc.pp.normalize_total(CAF, target_sum=1e4)
sc.pp.log1p(CAF)
sc.pp.highly_variable_genes(CAF, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(CAF, max_value=10)
sc.pp.neighbors(CAF, use_rep='X_scanvi_fix_linear',n_neighbors=35)
sc.pp.pca(CAF)
sc.tl.pca(CAF,n_comps=50)
sc.pp.neighbors(CAF, n_pcs=30)
sc.tl.umap(CAF,min_dist=0.5)
sc.pl.umap(CAF,color=['stromal_celltype'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/CAF.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/CAF.pdf',dpi = 300, bbox_inches='tight')
sc.tl.leiden(CAF,resolution=0.2)
sc.pl.umap(CAF,color=['leiden'],legend_loc='on data')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/CAF_leiden.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/CAF_leiden.pdf',dpi = 300, bbox_inches='tight')
# 提取所需列
df = CAF.obs[['leiden', 'preKNN_CytoTRACE2_Score']]

# 计算每个 celltype 组内的 preKNN_CytoTRACE2_Score 均值，并按均值排序
order = df.groupby('leiden')['preKNN_CytoTRACE2_Score'].median().sort_values(ascending=False).index
import seaborn as sns
# 设置配色
palette = sns.color_palette("Spectral", len(order))

plt.figure(figsize=(12, 6))

# 小提琴图（按 order 排序）
sns.violinplot(x='leiden', y='preKNN_CytoTRACE2_Score', data=df,
               order=order, palette=palette, inner=None, linewidth=1.2)

# 叠加箱型图（按 order 排序）
sns.boxplot(x='leiden', y='preKNN_CytoTRACE2_Score', data=df,
            order=order,
            showcaps=False,
            boxprops={'facecolor': 'none'},
            showfliers=False,
            whiskerprops={'linewidth': 1.5},
            width=0.4)

# 图像美化
plt.xticks(rotation=45, fontsize=10)
plt.title('preKNN_CytoTRACE2_Score by leiden (Ordered by Median)', fontsize=14)
plt.ylabel('preKNN_CytoTRACE2_Score', fontsize=12)
plt.xlabel('leiden', fontsize=12)
plt.tight_layout()

# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/stromal/'
plt.savefig(results_file + 'CAFviolin_boxplot_preKNN_CytoTRACE2_Score_sorted.png',
            dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAFviolin_boxplot_preKNN_CytoTRACE2_Score_sorted.pdf',
            dpi=300, bbox_inches='tight')
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 提取 GPNMB 表达量和分组信息
gene = 'GPNMB'

# 获取表达矩阵（假设使用的是 log1p 处理后的数据）
expr = stromal[:, gene].X.toarray().flatten()  # 如果是 dense matrix 可以直接用 .X

# 构建绘图用的 DataFrame
df = pd.DataFrame({
    'stromal_celltype': stromal.obs['stromal_celltype'].values,
    gene: expr
})

# 计算每个 celltype 组内的 median 表达量，并按从高到低排序
order = df.groupby('stromal_celltype')[gene].median().sort_values(ascending=False).index

# 设置配色方案
palette = sns.color_palette("Spectral", len(order))

# 设置画布大小
plt.figure(figsize=(12, 6))

# 小提琴图（按 order 排序）
sns.violinplot(x='stromal_celltype', y=gene, data=df,
               order=order, palette=palette, inner=None, linewidth=1.2)

# 叠加箱型图（按 order 排序）
sns.boxplot(x='stromal_celltype', y=gene, data=df,
            order=order,
            showcaps=False,
            boxprops={'facecolor': 'none'},
            showfliers=False,
            whiskerprops={'linewidth': 1.5},
            width=0.4)

# 图像美化
plt.xticks(rotation=45, fontsize=10)
plt.title(f'{gene} Expression by Stromal Celltype (Ordered by Median)', fontsize=14)
plt.ylabel(f'log1p({gene}) Expression', fontsize=12)
plt.xlabel('Stromal Celltype', fontsize=12)
plt.tight_layout()

# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/stromal/'
plt.savefig(results_file + f'violin_boxplot_{gene}_sorted_by_median.png',
            dpi=300, bbox_inches='tight')

sc.pl.umap(CAF,color=gene)
plt.savefig(results_file + f'CAFumap_{gene}.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + f'CAFumap_{gene}.pdf', dpi=300, bbox_inches='tight')
sc.tl.rank_genes_groups(CAF, method="t-test",groupby = 'leiden')
sc.pl.rank_genes_groups(CAF, n_genes=20, sharey=False)
plt.savefig(results_file + f'CAFrank_genes_groups_{gene}.png',
            dpi=300, bbox_inches='tight')
plt.savefig(results_file + f'CAFrank_genes_groups_{gene}.pdf',
            dpi=300, bbox_inches='tight')
sc.pl.umap(CAF,color='Tissue')
plt.savefig(results_file + 'CAF_umap_tissue.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAF_umap_tissue.pdf', dpi=300, bbox_inches='tight')
# 提供的 marker 基因列表
marker_genes = {
    "Human iCAF": [
        "ADAMTS4", "AGT", "APOE", "ARHGDI", "CCL19",
        "CCL21", "COLEC11", "CPE", "GEM", "GJA4", "GPX3",
        "HIGD1B", "IL6", "ISYNAI", "LHFP", "MAP1B",
        "MT1A", "NDUFA4L2", "PDK4", "RGS5"
    ],
    "Human myCAF": [
        "APOD", "CCL11", "COL1A1", "COL1A2", "COL3A1",
        "COL5A1", "COL6A3", "CTGF", "CTHRC1",
        "CYP1B1", "FN1", "INHBA", "ISLR", "LUM", "MMP14",
        "POSTN", "PTGDS", "SERPINFI", "SFRP2", "SPON2",
        "VCAN"
    ],
    "Human mesCAF": [
        "ANXA1", "ANXA2", "BDKRB1", "C19orf33", "C3",
        "CALB2", "CCDC80", "CFB", "CRABP2", "CXCL1",
        "CXCL6", "EFEMP1", "EGFL6", "EMP3", "EZR",
        "HMOX1", "HP", "HSPA6", "IFI27", "IGFBP6", "ITLN1",
        "KRT18", "KRT19", "KRT8", "LINC01133", "LOX",
        "MT1E", "MTIG", "MTIX", "MXRA5", "PDPN",
        "PLAZG2A", "PRG4", "PRSS23", "PTGIS", "RP11-572C15.6",
        "S100A10", "S100A16", "S100A6", "SAA1", "SAA2",
        "SERPINE2", "SH3BGR3", "SLC12A8", "SLPI", "TM4SF1"
    ],
    "Human vCAF": [
        "CCL8", "GJA4", "MHY11", "MCAM", "RGS5", "IL-6"
    ],
    "Human mCAF": [
        "COL5A1", "COL5A2", "COL6A3", "DCN", "FN1",
        "LUM", "POSTN", "VCAN"
    ],
    "Human apCAF": [
        "CD74", "HLA-DRA", "HLA-DRB1", "CCL21",
        "CXCL12"
    ],
    "Human eCAF": [
        "KRT19", "KRT8", "SAA1", "SLPI"
    ],
    "Human ICAF": [
        "APOA2", "FABP1", "FABP4", "FRZB", "GPX3"
    ],
    "Human iCAF": [
        "CXCL1", "C3", "C7", "FBLN1", "IGFI", "IGFBP6",
        "SAA1"
    ],
    "Human vCAF": [
        "ADIRF", "RGS5", "SPARCL1", "CRIP1",
        "NDUFA4L2", "MYH11", "MCAM", "PDK4", "FABP4",
        "TINAGL1"
    ],
    "Human mCAF": [
        "POSTN", "CTHRC1", "COL1A1", "COL6A3", "FN1",
        "COL3A1", "MMP14", "LUM", "SPON2", "COL5A1"
    ],
    "Human apCAF": [
        "IGFBP3", "CXCL12", "HLA-DRB1", "CD74", "HLA-DRA",
        "RBP1", "HLA-DPB1", "COLEC11"
    ],
    "Human imCAF": [
        "CCDC80", "C3", "SERPINFI", "PTGDS", "FBLN1",
        "IGF1", "CRABP2", "MMP2", "CTGF"
    ]
}

# 提取 leiden 分组信息
grouped_data = CAF.obs['leiden'].copy()

# 提取并排序所有 marker 基因（保留重复）
ordered_genes = []
gene_to_groups = {}  # 记录每个基因属于哪些 group
for group, genes in marker_genes.items():
    for gene in genes:
        if gene in CAF.var_names:
            ordered_genes.append(gene)
            if gene not in gene_to_groups:
                gene_to_groups[gene] = set()
            gene_to_groups[gene].add(group)

# 去重但保持顺序（用于提取表达矩阵）
unique_ordered_genes = []
seen = set()
for g in ordered_genes:
    if g not in seen:
        seen.add(g)
        unique_ordered_genes.append(g)

# 提取表达矩阵（log1p 标准化后）
expr_df = pd.DataFrame(CAF[:, unique_ordered_genes].X.toarray(),
                       index=CAF.obs_names,
                       columns=unique_ordered_genes)

# 按照 leiden 分组计算平均表达
heatmap_data = expr_df.groupby(grouped_data).mean().T  # Gene x Cluster

# 构建颜色条带（每个基因对应的颜色）
gene_colors = []
legend_elements = []
group_to_color = {}
color_palette = sns.color_palette("tab10", n_colors=len(marker_genes))

for group, color in zip(marker_genes.keys(), color_palette):
    group_to_color[group] = color

# 对每个基因，如果有多个 group，取第一个作为颜色显示
for gene in heatmap_data.index:
    groups = gene_to_groups.get(gene, [])
    if len(groups) > 0:
        first_group = list(groups)[0]
        gene_colors.append(group_to_color[first_group])
    else:
        gene_colors.append((0.8, 0.8, 0.8))  # 默认灰色

# 构建 legend
import matplotlib.patches as mpatches
for group, color in group_to_color.items():
    legend_elements.append(mpatches.Patch(color=color, label=group))

# 绘制 clustermap（带注释）
g = sns.clustermap(
    heatmap_data,
    cmap='viridis',
    z_score=0,  # 按基因标准化
    row_cluster=False,  # 不聚类基因
    col_cluster=False,  # 不聚类样本
    yticklabels=True,
    xticklabels=True,
    figsize=(10, 20),
    row_colors=[gene_colors],
    cbar_kws={'label': 'Z-score Expression'}
)

# 添加图例
g.ax_row_dendrogram.set_visible(False)
g.fig.subplots_adjust(right=0.7)
g.fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.2, 0.8))

# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/stromal/'
plt.savefig(results_file + 'CAFclustermap_with_gene_annotations.png', dpi=300, bbox_inches='tight')
sc.pl.umap(CAF,color = 'Origin')
plt.savefig(results_file + 'CAFOrigin.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAFOrigin.pdf', dpi=300, bbox_inches='tight')
CAF.obs['leiden'].cat.categories
sc.tl.rank_genes_groups(CAF, 'leiden', method='t-test')

sc.pl.rank_genes_groups(CAF, n_genes=25, sharey=False)
plt.savefig(results_file + 'CAFrank_genes_groups.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAFrank_genes_groups.pdf', dpi=300, bbox_inches='tight')

sc.pl.umap(CAF,color = 'ATP1B3')
plt.savefig(results_file + 'CAFumap_ATP1B3.png', dpi=300, bbox_inches='tight')
map = {'0':'POSTN+ mCAF',
       '1':'Alveolar CAF',
       '2':'COL15A1+ mCAF',
       '3':'MSLN+ iCAF',
       '4':'PI16+ proCAF',
       '5':'LAMA2+ myCAF',
       '6':'COL11A1+ myCAF'}
CAF.obs['CAF_type'] = CAF.obs['leiden'].map(map).astype('category')
marker_genes = {'POSTN+ mCAF':['POSTN',"COL1A1", "COL1A2", "COL3A1"],
                'Alveolar CAF' : ['A2M', 'LIMCH1', 'SCN7A', 'MYH10'],
                'COL15A1+ mCAF':['COL5A2', 'COL3A1', 'THBS2', 'COL15A1'],
                'MSLN+ iCAF':['MSLN','VEGFA','AKR1C1','AKR1C2'],
                'PI16+ proCAF':['PI16','IGFBP6','CST3','FBLN1'],
                'LAMA2+ myCAF':['LAMA2','APOD','NTRK3','CFH'],
                'COL11A1+ myCAF':['COL11A1','FOXP1','FOXP2','LEPR']}
# 提取 leiden 分组信息
grouped_data = CAF.obs['CAF_type'].copy()

# 提取并排序所有 marker 基因（保留重复）
ordered_genes = []
gene_to_groups = {}  # 记录每个基因属于哪些 group
for group, genes in marker_genes.items():
    for gene in genes:
        if gene in CAF.var_names:
            ordered_genes.append(gene)
            if gene not in gene_to_groups:
                gene_to_groups[gene] = set()
            gene_to_groups[gene].add(group)

# 去重但保持顺序（用于提取表达矩阵）
unique_ordered_genes = []
seen = set()
for g in ordered_genes:
    if g not in seen:
        seen.add(g)
        unique_ordered_genes.append(g)

# 提取表达矩阵（log1p 标准化后）

expr_df = pd.DataFrame(CAF[:, unique_ordered_genes].X.toarray(),
                       index=CAF.obs_names,
                       columns=unique_ordered_genes)

# 按照 leiden 分组计算平均表达
heatmap_data = expr_df.groupby(grouped_data).mean().T  # Gene x Cluster

# 构建颜色条带（每个基因对应的颜色）
gene_colors = []
legend_elements = []
group_to_color = {}
color_palette = sns.color_palette("tab10", n_colors=len(marker_genes))

for group, color in zip(marker_genes.keys(), color_palette):
    group_to_color[group] = color

# 对每个基因，如果有多个 group，取第一个作为颜色显示
for gene in heatmap_data.index:
    groups = gene_to_groups.get(gene, [])
    if len(groups) > 0:
        first_group = list(groups)[0]
        gene_colors.append(group_to_color[first_group])
    else:
        gene_colors.append((0.8, 0.8, 0.8))  # 默认灰色

# 构建 legend
import matplotlib.patches as mpatches
for group, color in group_to_color.items():
    legend_elements.append(mpatches.Patch(color=color, label=group))

# 绘制 clustermap（带注释）
g = sns.clustermap(
    heatmap_data,
    cmap='viridis',
    z_score=0,  # 按基因标准化
    row_cluster=False,  # 不聚类基因
    col_cluster=False,  # 不聚类样本
    yticklabels=True,
    xticklabels=True,
    figsize=(10, 15),
    row_colors=[gene_colors],
    cbar_kws={'label': 'Z-score Expression'}
)

# 添加图例
g.ax_row_dendrogram.set_visible(False)
g.fig.subplots_adjust(right=0.7)
g.fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.2, 0.8))

# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/stromal/'
plt.savefig(results_file + 'CAFclustermap_with_gene_annotations.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAFclustermap_with_gene_annotations.pdf', dpi=300, bbox_inches='tight')
# 假设你已经有了 CAF 对象，并且已经定义了 CAF_type
df = CAF.obs[['Tissue', 'CAF_type']].copy()

# 固定 Tissue 的顺序
desired_order = ['tumor_middle', 'tumor_edge', 'normal_adjacent', 'tumor_primary', 'tumor_metastasis']
df['Tissue'] = pd.Categorical(df['Tissue'], categories=desired_order, ordered=True)

# 构建交叉表并计算比例
count_df = pd.crosstab(index=df['Tissue'], columns=df['CAF_type'])
prop_df = count_df.div(count_df.sum(axis=1), axis=0)

# 设置画布
plt.figure(figsize=(12, 6))

# 使用固定的配色方案（根据 CAF_type 排序）
colors = sns.color_palette("Spectral", n_colors=len(prop_df.columns))
ax = prop_df.plot(kind='bar', stacked=True, color=colors, edgecolor='black')

# 添加百分比标签（可选）
for i, (idx, row) in enumerate(prop_df.iterrows()):
    height = 0
    for col_name, val in row.items():
        if val > 0.05:  # 只对大于 5% 的部分添加文字
            ax.text(i, height + val / 2., f'{val:.1%}', ha='center', va='center', fontsize=9)
            height += val

# 美化图表
plt.title('Proportion of CAF Subtypes by Tissue (Specified Order)', fontsize=14)
plt.xlabel('Tissue')
plt.ylabel('Proportion')
plt.xticks(rotation=45)
plt.legend(title='CAF Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/stromal/'
plt.savefig(results_file + 'CAF_tissue_caf_type_stacked_barplot_ordered.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAF_tissue_caf_type_stacked_barplot_ordered.pdf', dpi=300, bbox_inches='tight')
import math

# 获取所有 tissue 类别（已排序）
tissues = desired_order

# 设置子图行列数
cols = 3
rows = math.ceil(len(tissues) / cols)

# 创建画布
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

# 绘图
for i, tissue in enumerate(tissues):
    subset = df[df['Tissue'] == tissue]
    counts = subset['CAF_type'].value_counts()
    
    counts.plot(kind='pie', ax=axes[i], autopct='%1.1f%%',
                title=f'{tissue}', legend=False,
                colors=sns.color_palette("Spectral", len(counts)))
    axes[i].set_ylabel('')

# 隐藏多余子图
for j in range(i+1, rows*cols):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig(results_file + 'CAF_piecharts_by_Tissue_ordered.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAF_piecharts_by_Tissue_ordered.pdf', dpi=300, bbox_inches='tight')
# 提取所需列
df = CAF.obs[['CAF_type', 'preKNN_CytoTRACE2_Score']]

# 计算每个 celltype 组内的 preKNN_CytoTRACE2_Score 均值，并按均值排序
order = df.groupby('CAF_type')['preKNN_CytoTRACE2_Score'].median().sort_values(ascending=False).index
import seaborn as sns
# 设置配色
palette = sns.color_palette("Spectral", len(order))

plt.figure(figsize=(12, 6))

# 小提琴图（按 order 排序）
sns.violinplot(x='CAF_type', y='preKNN_CytoTRACE2_Score', data=df,
               order=order, palette=palette, inner=None, linewidth=1.2)

# 叠加箱型图（按 order 排序）
sns.boxplot(x='leiden', y='preKNN_CytoTRACE2_Score', data=df,
            order=order,
            showcaps=False,
            boxprops={'facecolor': 'none'},
            showfliers=False,
            whiskerprops={'linewidth': 1.5},
            width=0.4)

# 图像美化
plt.xticks(rotation=45, fontsize=10)
plt.title('preKNN_CytoTRACE2_Score by CAF_type (Ordered by Median)', fontsize=14)
plt.ylabel('preKNN_CytoTRACE2_Score', fontsize=12)
plt.xlabel('CAF_type', fontsize=12)
plt.tight_layout()

# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/stromal/'
plt.savefig(results_file + 'CAF_typeviolin_boxplot_preKNN_CytoTRACE2_Score_sorted.png',
            dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAF_typeviolin_boxplot_preKNN_CytoTRACE2_Score_sorted.pdf',
            dpi=300, bbox_inches='tight')
sc.tl.diffmap(CAF)
sc.pl.diffmap(CAF,color = 'CAF_type')
plt.savefig(results_file + 'CAF_type_diffmap.png', dpi=300, bbox_inches='tight')

CAF.obsm['X_diffmap_'] = CAF.obsm['X_diffmap'][:,1:]
sc.pl.embedding(CAF,'diffmap_',color=['CAF_type'])
plt.savefig(results_file + 'CAF_type_diffmap.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAF_type_diffmap.pdf', dpi=300, bbox_inches='tight')
sc.pp.neighbors(CAF, n_neighbors=20, use_rep='X_diffmap_')
sc.tl.draw_graph(CAF)
sc.pl.draw_graph(CAF, color='CAF_type')
plt.savefig(results_file + 'CAF_type_draw_graph.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAF_type_draw_graph.pdf', dpi=300, bbox_inches='tight')
sc.tl.paga(CAF, groups='CAF_type')
sc.pl.paga(CAF, color='CAF_type')
plt.savefig(results_file + 'CAF_paga_CAF_type.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAF_paga_CAF_type.pdf', dpi=300, bbox_inches='tight')
CAF.uns["iroot"] = np.flatnonzero(CAF.obs["CAF_type"] == "PI16+ proCAF")[0]

sc.tl.dpt(CAF)
sc.pl.draw_graph(CAF, color=["dpt_pseudotime"], legend_loc="on data")
plt.savefig(results_file + 'CAF_dpt_pseudotime.png',dpi = 300, bbox_inches='tight')
plt.savefig(results_file + 'CAF_dpt_pseudotime.pdf',dpi = 300, bbox_inches='tight')
sc.pl.paga(CAF,color=["dpt_pseudotime"])
plt.savefig(results_file + 'CAF_dpt_pseudotime_paga.png',dpi = 300, bbox_inches='tight')
plt.savefig(results_file + 'CAF_dpt_pseudotime_paga.pdf',dpi = 300, bbox_inches='tight')
# 计算每个 leiden group 的 pseudotime 平均表达
# 添加伪时间排序
CAF.obs = CAF.obs.sort_values(by='dpt_pseudotime')

# 构建一个 DataFrame，包含 pseudotime 和 所有要绘图的基因 + leiden 分组
gene_df = CAF[:, unique_ordered_genes].to_df()
plot_df = gene_df.copy()
plot_df['dpt_pseudotime'] = CAF.obs['dpt_pseudotime'].values
plot_df['CAF_type'] = CAF.obs['CAF_type'].astype(str).values

# 设置画布大小和颜色
n_clusters = len(plot_df['CAF_type'].unique())
colors = sns.color_palette("Spectral", n_clusters)

ordered_df = []



# 每个基因单独画一张图，或在一个图中叠加多个基因
for gene in unique_ordered_genes:
    plt.figure(figsize=(10, 6))
    
    for i, cluster in enumerate(sorted(plot_df['CAF_type'].unique())):
        sub_df = plot_df[plot_df['CAF_type'] == cluster]
        sns.lineplot(data=sub_df, x='dpt_pseudotime', y=gene, color=colors[i], label=f'CAF_type {cluster}')
    
    # 美化设置
    plt.title(f'{gene} Expression Along Pseudotime (Colored by CAF_type Cluster)')
    plt.xlabel('Pseudotime')
    plt.ylabel('Expression (log1p)')
    plt.legend(title='CAF_type Cluster')
    plt.tight_layout()

    # 保存图像
    results_file = '/home/data/sdzl14/NSCLC/zong/fig/stromal/'
    plt.savefig(results_file + f'CAFpseudotime_{gene}_by_CAF_type.png', dpi=300, bbox_inches='tight')
    plt.close()
CAF.obs['Timepoint'].value_counts()
CAF.write_h5ad('/home/data/sdzl14/NSCLC/zong/CAF.h5ad')
# 假设你已经有了 CAF 对象，并且已经定义了 CAF_type
df = CAF.obs[['Timepoint', 'CAF_type']].copy()

# 固定 Tissue 的顺序
desired_order = ['Pre', 'Post']
df['Timepoint'] = pd.Categorical(df['Timepoint'], categories=desired_order, ordered=True)

# 构建交叉表并计算比例
count_df = pd.crosstab(index=df['Timepoint'], columns=df['CAF_type'])
prop_df = count_df.div(count_df.sum(axis=1), axis=0)

# 设置画布
plt.figure(figsize=(12, 6))

# 使用固定的配色方案（根据 CAF_type 排序）
colors = sns.color_palette("Spectral", n_colors=len(prop_df.columns))
ax = prop_df.plot(kind='bar', stacked=True, color=colors, edgecolor='black')

# 添加百分比标签（可选）
for i, (idx, row) in enumerate(prop_df.iterrows()):
    height = 0
    for col_name, val in row.items():
        if val > 0.05:  # 只对大于 5% 的部分添加文字
            ax.text(i, height + val / 2., f'{val:.1%}', ha='center', va='center', fontsize=9)
            height += val

# 美化图表
plt.title('Proportion of CAF Subtypes by Timepoint (Specified Order)', fontsize=14)
plt.xlabel('Timepoint')
plt.ylabel('Proportion')
plt.xticks(rotation=45)
plt.legend(title='CAF Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/stromal/'
plt.savefig(results_file + 'CAF_Timepoint_caf_type_stacked_barplot_ordered.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAF_Timepoint_caf_type_stacked_barplot_ordered.pdf', dpi=300, bbox_inches='tight')
import math

# 获取所有 Timepoint 类别（已排序）
Timepoint = desired_order

# 设置子图行列数
cols = 3
rows = math.ceil(len(Timepoint) / cols)

# 创建画布
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

# 绘图
for i, Timepoint in enumerate(Timepoint):
    subset = df[df['Timepoint'] == Timepoint]
    counts = subset['CAF_type'].value_counts()
    
    counts.plot(kind='pie', ax=axes[i], autopct='%1.1f%%',
                title=f'{Timepoint}', legend=False,
                colors=sns.color_palette("Spectral", len(counts)))
    axes[i].set_ylabel('')

# 隐藏多余子图
for j in range(i+1, rows*cols):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig(results_file + 'CAF_piecharts_by_Timepoint_ordered.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'CAF_piecharts_by_Timepoint_ordered.pdf', dpi=300, bbox_inches='tight')
# 计算差异基因
sc.tl.rank_genes_groups(CAF, groupby='CAF_type', method='t-test')
# 获取 rank_genes_groups 结果
results = CAF.uns['rank_genes_groups']
# 获取 groupby 字段名称（即 CAF_type 中的所有类别）
groups = list(results['names'].dtype.fields.keys())

# 构建一个 DataFrame，每列是不同 group 的 top 基因
max_gene_count = len(results['names'][groups[0]])  # 假设每组基因数一致

data = {}
for group in groups:
    data[group] = results['names'][group][:]  # 提取该组所有基因名

# 获取 rank_genes_groups 的结果并转换为 DataFrame
df_rank_genes = sc.get.rank_genes_groups_df(CAF,log2fc_min=1,group=None,pval_cutoff=0.05)


# 提取字段：genes、logfoldchanges、pvals、pvals_adj、scores
df = pd.DataFrame({
    'gene': results['names'][group],
    'logfoldchanges': results['logfoldchanges'][group],
    'pvals': results['pvals'][group],
    'pvals_adj': results['pvals_adj'][group],
    'scores': results['scores'][group]
})
df_rank_genes.to_csv(results_file + 'rank_genes_list.csv', index=False)

# 获取每个 cluster 的 top 基因
result = CAF.uns['rank_genes_groups']
groups = result['names'].dtype.names
n_top_genes = 2000  # 可根据需要调整

# 提取每个 cluster 的 top 基因
marker_genes_dict = {
    group: [result['names'][group][i] for i in range(n_top_genes)]
    for group in groups
}
import gseapy as gp
# 选择一个 cluster 的 top genes 举例
genes_of_interest = marker_genes_dict['POSTN+ mCAF']  # 替换为你感兴趣的 cluster 名称
import pandas as pd
immune = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/immune.subset.h5ad')
immune = immune.copy()
immune.obs['immune_celltype_coarse'].cat.categories
T  = ['CD8','CD4','Cycling T','Treg']
DC = ['cDC1','cDC2','mregDC','pDC']
T = immune[immune.obs['immune_celltype_coarse'].isin(T)]
DC = immune[immune.obs['immune_celltype_coarse'].isin(DC)]
B = ['Cycling B','GC B','Memory B','Naive B','Plasma']
B = immune[immune.obs['immune_celltype_coarse'].isin(B)]
Macro = immune[immune.obs['immune_celltype_coarse'] == 'Macro']
# 假设你已经有了 CAF 对象，并且已经定义了 CAF_type
df = Macro.obs[['Pathtype', 'immune_celltype']].copy()

# 固定 Tissue 的顺序
desired_order = ['LUSC','LUAD','NSCLC']

df['Pathtype'] = pd.Categorical(df['Pathtype'], categories=desired_order, ordered=True)

# 构建交叉表并计算比例
count_df = pd.crosstab(index=df['Pathtype'], columns=df['immune_celltype'])
prop_df = count_df.div(count_df.sum(axis=1), axis=0)

# 设置画布
plt.figure(figsize=(12, 6))

# 使用固定的配色方案（根据 immune_celltype 排序）
colors = sns.color_palette("Spectral", n_colors=len(prop_df.columns))
ax = prop_df.plot(kind='bar', stacked=True, color=colors, edgecolor='black')

# 添加百分比标签（可选）
for i, (idx, row) in enumerate(prop_df.iterrows()):
    height = 0
    for col_name, val in row.items():
        if val > 0.05:  # 只对大于 5% 的部分添加文字
            ax.text(i, height + val / 2., f'{val:.1%}', ha='center', va='center', fontsize=9)
            height += val

# 美化图表
plt.title('Proportion of macro Subtypes by Pathtype (Specified Order)', fontsize=14)
plt.xlabel('Tissue')
plt.ylabel('Proportion')
plt.xticks(rotation=45)
plt.legend(title='Macro Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/immune/'
plt.savefig(results_file + 'macro_Pathtype_macro_type_stacked_barplot_ordered.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'macro_Pathtype_macro_type_stacked_barplot_ordered.pdf', dpi=300, bbox_inches='tight')
import math

# 获取所有 tissue 类别（已排序）
tissues = desired_order

# 设置子图行列数
cols = 3
rows = math.ceil(len(tissues) / cols)

# 创建画布
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

# 绘图
for i, tissue in enumerate(tissues):
    subset = df[df['Tissue'] == tissue]
    counts = subset['immune_celltype'].value_counts()
    
    counts.plot(kind='pie', ax=axes[i], autopct='%1.1f%%',
                title=f'{tissue}', legend=False,
                colors=sns.color_palette("Spectral", len(counts)))
    axes[i].set_ylabel('')

# 隐藏多余子图
for j in range(i+1, rows*cols):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig(results_file + 'macro_piecharts_by_Tissue_ordered.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'macro_piecharts_by_Tissue_ordered.pdf', dpi=300, bbox_inches='tight')
# 假设你已经有了 CAF 对象，并且已经定义了 CAF_type
df = T.obs[['Tissue', 'immune_celltype']].copy()

# 固定 Tissue 的顺序
desired_order = ['tumor_middle', 'tumor_edge', 'normal_adjacent', 'tumor_primary', 'tumor_metastasis']
df['Tissue'] = pd.Categorical(df['Tissue'], categories=desired_order, ordered=True)

# 构建交叉表并计算比例
count_df = pd.crosstab(index=df['Tissue'], columns=df['immune_celltype'])
prop_df = count_df.div(count_df.sum(axis=1), axis=0)

# 设置画布
plt.figure(figsize=(12, 6))

# 使用固定的配色方案（根据 immune_celltype 排序）
colors = sns.color_palette("Spectral", n_colors=len(prop_df.columns))
ax = prop_df.plot(kind='bar', stacked=True, color=colors, edgecolor='black')

# 添加百分比标签（可选）
for i, (idx, row) in enumerate(prop_df.iterrows()):
    height = 0
    for col_name, val in row.items():
        if val > 0.05:  # 只对大于 5% 的部分添加文字
            ax.text(i, height + val / 2., f'{val:.1%}', ha='center', va='center', fontsize=9)
            height += val

# 美化图表
plt.title('Proportion of T Subtypes by Tissue (Specified Order)', fontsize=14)
plt.xlabel('Tissue')
plt.ylabel('Proportion')
plt.xticks(rotation=45)
plt.legend(title='T Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/immune/'
plt.savefig(results_file + 'T_tissue_T_type_stacked_barplot_ordered.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'T_tissue_T_type_stacked_barplot_ordered.pdf', dpi=300, bbox_inches='tight')
import math

# 获取所有 tissue 类别（已排序）
tissues = desired_order

# 设置子图行列数
cols = 3
rows = math.ceil(len(tissues) / cols)

# 创建画布
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

# 绘图
for i, tissue in enumerate(tissues):
    subset = df[df['Tissue'] == tissue]
    counts = subset['immune_celltype'].value_counts()
    
    counts.plot(kind='pie', ax=axes[i], autopct='%1.1f%%',
                title=f'{tissue}', legend=False,
                colors=sns.color_palette("Spectral", len(counts)))
    axes[i].set_ylabel('')

# 隐藏多余子图
for j in range(i+1, rows*cols):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig(results_file + 'T_piecharts_by_Tissue_ordered.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'Tpiecharts_by_Tissue_ordered.pdf', dpi=300, bbox_inches='tight')
# 假设你已经有了 CAF 对象，并且已经定义了 CAF_type

df = immune.obs[['Timepoint', 'immune_celltype_coarse']].copy()
B.obs['Pathtype']
# 固定 Tissue 的顺序
desired_order =['LUAD', 'LUSC', 'NSCLC', 'NSCLC NOS']
df['Timepoint'] = pd.Categorical(df['Timepoint'], categories=desired_order, ordered=True)

# 构建交叉表并计算比例
count_df = pd.crosstab(index=df['Timepoint'], columns=df['immune_celltype_coarse'])
prop_df = count_df.div(count_df.sum(axis=1), axis=0)
#创建新的 DataFrame 来存储处理后的比例数据
prop_df_with_other = pd.DataFrame(index=prop_df.index, columns=prop_df.columns.tolist() + ['Other'])

# 遍历每一行（每个 Timepoint "Other"
threshold = 0.05
for tissue in prop_df.index:
    row = prop_df.loc[tissue]
    # 找出比例大于等于阈值的细胞类型
    valid_types = row[row >= threshold]
    # 计算 Other 的比例
    other_ratio = 1 - valid_types.sum()
    # 将有效类型和 Other 拼接为新的行数据
    new_row = pd.concat([valid_types, pd.Series(other_ratio, index=['Other'])])
    prop_df_with_other.loc[tissue] = new_row

# 填充缺失值（当所有类型都低于阈值时可能会出现 NaN）
prop_df_with_other.fillna(0, inplace=True)
# 删除所有值为 0 的列
prop_df_with_other = prop_df_with_other.loc[:, (prop_df_with_other != 0).any(axis=0)]
# 现在可以使用 prop_df_with_other 继续绘图
# 设置画布
plt.figure(figsize=(12, 12))

# 使用固定的配色方案（根据 prop_df_with_other 列的数量）
colors = sns.color_palette("Spectral", n_colors=len(prop_df_with_other.columns))
ax = prop_df_with_other.plot(kind='bar', stacked=True, color=colors, edgecolor='black')

# 添加百分比标签
for i, (idx, row) in enumerate(prop_df_with_other.iterrows()):
    height = 0
    for col_name, val in row.items():
        if val > 0.05:  # 对大于 5% 的部分添加文字
            ax.text(i, height + val / 2., f'{val:.1%}', ha='center', va='center', fontsize=9)
            height += val

# 美化图表
plt.title('Proportion of B Subtypes by Timepoint (Specified Order)', fontsize=14)
plt.xlabel('Timepoint')
plt.ylabel('Proportion')
plt.xticks(rotation=45)
plt.legend(title='Immune Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()



# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/immune/'
plt.savefig(results_file + 'Immune_Timepoint_immune_type_stacked_barplot_ordered.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'Immune_Timepoint_immune_type_stacked_barplot_ordered.pdf', dpi=300, bbox_inches='tight')
# 获取所有 Tissue 类别（已排序）
Timepoint = desired_order

# 设置子图行列数
cols = 3
rows = math.ceil(len(Timepoint) / cols)

# 创建画布
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

# 统一配色方案（基于所有可能的 celltypes）
all_celltypes = prop_df_with_other.columns.tolist()
colors = sns.color_palette("Spectral", len(all_celltypes))

# 绘图
for i, tissue in enumerate(Timepoint):
    subset = prop_df_with_other.loc[[tissue]].T.dropna()
    subset = subset[subset[tissue] > 0]  # 过滤掉 0 值
    labels = subset.index.tolist()
    values = subset[tissue].tolist()
    
    # 绘制饼图
    axes[i].pie(
        values,
        labels=labels,
        autopct='%1.1f%%',
        colors=[colors[all_celltypes.index(label)] for label in labels],
        textprops={'fontsize': 10}
    )
    axes[i].set_title(tissue)
    axes[i].axis('equal')  # 确保饼图为圆形

# 隐藏多余子图
for j in range(i + 1, rows * cols):
    fig.delaxes(axes[j])

plt.tight_layout()
results_file = '/home/data/sdzl14/NSCLC/zong/fig/immune/'
plt.savefig(results_file + 'B_piecharts_by_Timepoint_ordered_with_Other.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'B_piecharts_by_Timepoint_ordered_with_Other.pdf', dpi=300, bbox_inches='tight')
sc.pl.umap(immune ,color=['immune_celltype_coarse'])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/immune_celltype_coarse.png',dpi = 300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/immune/immune_celltype_coarse.pdf',dpi = 300, bbox_inches='tight')
