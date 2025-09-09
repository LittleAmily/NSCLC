from natsort import order_by_index
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
malignant = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/malignant.h5ad')
malignant = malignant.copy()
LUSC = malignant[malignant.obs['Pathtype'] == 'LUSC']
sc.tl.rank_genes_groups(LUSC, 'Tissue', method='t-test')
sc.pl.rank_genes_groups(LUSC, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_intergrated_rank_genes_groups.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_intergrated_rank_genes_groups.png",dpi = 300, bbox_inches='tight')
LUSC.obs['Tissue'].value_counts()
map = {'effusion':'effusion', 
       'normal_adjacent':'normal_adjacent', 
       'normal_adjancent':'normal_adjacent',
       'tumor_edge':'tumor_edge', 
       'tumor_metastasis':'tumor_metastasis', 
       'tumor_metastatis':'tumor_metastasis', 
       'tumor_middle':'tumor_middle',
       'tumor_primary':'tumor_primary'}
LUSC.obs['Tissue'] = (
    LUSC.obs['Tissue']
    .map(map)
    .astype('category')
)
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
order = ['tumor_middle','tumor_edge', 'normal_adjacent', 'tumor_metastasis', 'effusion']
LUSC_ = LUSC[LUSC.obs['Tissue'].isin(order)]
sc.pp.neighbors(LUSC_,use_rep='X_scANVI')
sc.tl.umap(LUSC_)
sc.tl.leiden(LUSC_, resolution=0.2)
sc.pl.umap(LUSC_, color="leiden")
LUSC_.obs['leiden'].cat.categories
map = {'0':'0', 
       '1':'1', 
       '2':'2', 
       '3':'3', 
       '4':'4', 
       '5':'5', 
       '6':'6', 
       '7':'7', 
       '8':'5'}
LUSC_.obs['leiden'] = (
    LUSC_.obs['leiden']
    .map(map)
    .astype('category')
)   
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap.pdf",dpi = 300, bbox_inches='tight')
sc.pl.umap(LUSC_, color="Tissue")
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Tissue.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Tissue.pdf",dpi = 300, bbox_inches='tight')
sc.pl.umap(LUSC_, color="Patient")
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Patient.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Patient.pdf",dpi = 300, bbox_inches='tight')
sc.pl.umap(LUSC_, color="Origin")
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Origin.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Origin.pdf",dpi = 300, bbox_inches='tight')
sc.pl.umap(LUSC_, color="Timepoint")
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Timepoint.png",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_umap_Timepoint.pdf",dpi = 300, bbox_inches='tight')
sc.tl.rank_genes_groups(LUSC_, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(LUSC_, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_intergrated_rank_genes_groups_.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_intergrated_rank_genes_groups_.png",dpi = 300, bbox_inches='tight')
sc.tl.rank_genes_groups(LUSC_, 'Tissue', method='logreg')
sc.pl.rank_genes_groups(LUSC_, n_genes=25, sharey=False)
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_rank_genes_groups_.pdf",dpi = 300, bbox_inches='tight')
plt.savefig("/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_rank_genes_groups_.png",dpi = 300, bbox_inches='tight')
# 假设你的adata是AnnData对象
import pandas as pd
from anndata import AnnData
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

df_result = calculate_cell_percentage(LUSC_, leiden_col='leiden', output_path='/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_Origincell_counts_percentage.csv')
# 验证每个组织的百分比总和应为100%
check = df_result.groupby('Tissue')['percentage'].sum()
print("\n组织类型百分比验证:")
print(check.round(2))
LUSC_.obs['Dataset'].cat.categories
LUSC = ['Goveia_Carmeliet_2020', 'He_Fan_2021', 'Kim_Lee_2020',
       'Lambrechts_Thienpont_2018_6149v2', 'Lambrechts_Thienpont_2018_6653',
       'Laughney_Massague_2020', 'Maynard_Bivona_2020', 'Peilin_Wang_2025', 'UKIM-V']
LUSC_ = LUSC_[LUSC_.obs['Dataset'].isin(LUSC)]
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
    tissue_order = ['tumor_middle','tumor_edge', 'normal_adjacent', 'tumor_metastasis']
    
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
# 设置目标基因和组织类型
target_genes = ['YBX1','PTPRC','PPDPF','FTL']
tissue_order = ['tumor_middle','tumor_edge', 'normal_adjacent', 'tumor_metastasis']
LUSC_ = LUSC
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
gene_palette = ['#f9d580','#85a2d2','#e3716e','#eca680','#7ac7e2','#f7df87','#54beaa','#2983b1','#882E71'] 

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
plt.title("Gene Expression Across Tissue Types", 
         fontsize=14, pad=20)
plt.xlabel("Tissue Category", fontsize=12)
plt.ylabel("Log Normalized Expression", fontsize=12)
plt.xticks(rotation=45, ha='right', fontsize=11)
plt.legend(title='Gene', frameon=False, bbox_to_anchor=(1, 0.5))

# 保存输出
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_expression-.pdf', 
           bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_expression-.png', 
           dpi=300, bbox_inches='tight')
plt.close()
target_genes = ['SFTPC','RPL17','B2M','RPL41','RPS6','NME2','EEF1A1','APOC1','SCGB1A1','TUBA1B','PTMA']
tissue_order = ['tumor_middle','tumor_edge', 'normal_adjacent', 'tumor_metastasis']

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
gene_palette =['#fcf1f0', '#fccccb','#bdb5e1','#b0d992','#f9d580','#85a2d2','#e3716e','#eca680','#7ac7e2','#f7df87','#54beaa','#2983b1','#882E71']  # 颜色数量需与Tissue类别数一致


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
plt.title("S100A Gene Expression Across Tissue Types", 
         fontsize=14, pad=20)
plt.xlabel("Tissue Category", fontsize=12)
plt.ylabel("Log Normalized Expression", fontsize=12)
plt.xticks(rotation=45, ha='right', fontsize=11)
plt.legend(title='Gene', frameon=False, bbox_to_anchor=(1, 0.5))

# 保存输出
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_expression__.pdf', 
           bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_expression__.png', 
           dpi=300, bbox_inches='tight')
LUSC_.var['mt'] = LUSC_.var_names.str.startswith('MT-')
# ribosomal genes
LUSC_.var['ribo'] = LUSC_.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
LUSC_.var['hb'] = LUSC_.var_names.str.contains(('HBA1|HBA2|HBB|HBD|HBE1|HBG1|HBG2|HBM|HBQ1|HBZ'))

# 计算质控指标
sc.pp.calculate_qc_metrics(LUSC_, qc_vars=['mt', 'ribo', 'hb'],layer='counts',
                           percent_top=None, log1p=False, inplace=True)
sc.pl.violin(LUSC_, ['pct_counts_mt', 'pct_counts_hb', 'pct_counts_ribo'], jitter=0.4,
 groupby = 'Dataset', rotation= 90)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_qc.pdf', 
           bbox_inches='tight')
sc.pl.violin(LUSC_, ['n_genes_by_counts','total_counts'], groupby = 'leiden', rotation= 90)
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/LUSC_leiden_qc2.pdf', 
           bbox_inches='tight')
LUSC_.write_h5ad("/home/data/sdzl14/NSCLC/zong/LUSC_qc.h5ad")
LUSC = sc.read_h5ad("/home/data/sdzl14/NSCLC/zong/LUSC.h5ad")
malignant = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/malignant_integrated.h5ad')