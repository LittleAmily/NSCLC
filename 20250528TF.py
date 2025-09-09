import pandas as pd
import anndata
import seaborn as sns
import matplotlib.pyplot as plt

# 步骤1: 读取Excel文件
hgnc_df = pd.read_excel('/home/data/sdzl14/NSCLC/zong/TranscriptionFactor.xlsx')
transcription_factors = hgnc_df['HGNC symbol'].tolist()

# 步骤2: 加载AnnData对象
adata = anndata.read_h5ad('/home/data/sdzl14/NSCLC/zong/malignant_integrated.h5ad')
adata = adata.copy()
adata.obs['Tissue'].cat.categories
# 步骤3: 过滤出转录因子相关的基因
tf_mask = adata.var_names.isin(transcription_factors)
adata_tf = adata[:, tf_mask]

# 步骤4: 按leiden分组并计算平均表达量
Tissue_avg_exp = adata_tf.to_df().groupby(adata.obs['Tissue']).mean()

# 步骤5: 获取每个leiden中Top5转录因子
top_tfs_per_Tissue = {}
for Tissue in Tissue_avg_exp.index:
    top_tfs = Tissue_avg_exp.loc[Tissue].sort_values(ascending=False).head(5)
    top_tfs_per_Tissue[Tissue] = top_tfs.index.tolist()

# 步骤6: 可视化不同leiden中的Top5转录因子
# 指定 leiden 顺序并过滤
custom_Tissue_order = ['tumor_middle','tumor_edge','normal_adjacent','tumor_metastasis']
Tissue_avg_exp_filtered = Tissue_avg_exp.reindex(custom_Tissue_order)

# 获取每个 Tissue 的 Top5 TF（各自独立）
top5_tfs_per_Tissue = {
    Tissue: Tissue_avg_exp_filtered.loc[Tissue].sort_values(ascending=False).head(5).index.tolist()
    for Tissue in custom_Tissue_order
}

# 合并所有 Tissue 的 Top5 TF 并去重
combined_top_tfs = list(set(tf for tfs in top5_tfs_per_Tissue.values() for tf in tfs))

# 构建热图数据：行 = Tissue，列 = combined_top_tfs
selected_tfs = ["E2F1", "E2F2", "E2F3","GLI1", "GLI2", "GLI3","FOS","FOSL1","JUN","RBPJ","MYC", "HES1", "HEY1"]
heatmap_data = Tissue_avg_exp_filtered[selected_tfs]

# 绘制热图
plt.figure(figsize=(6, 6))
sns.heatmap(heatmap_data, annot=True, cmap='viridis', linewidths=.5, cbar_kws={'label': 'Expression'})
plt.title('Expression of Combined Top Transcription Factors Across Tissues')
plt.xlabel('Transcription Factor')
plt.ylabel('Tissue')
plt.xticks(rotation=45)
# 🔍 使用 seaborn.clustermap 实现聚类热图
g = sns.clustermap(
    data=heatmap_data,
    cmap='vlag',                   # 更温和的发散型配色（适合表达量数据）
    row_cluster=False,              # 对 Tissue（行）聚类
    col_cluster=False,              # 对 TF（列）聚类
    annot=True,                    # 显示数值
    fmt=".2f",                     # 注释保留两位小数
    linewidths=.5,                 # 格子间隔线
    annot_kws={"size": 14},        # 注释字体大小
    cbar_pos=(0.02, 0.8, 0.03, 0.15),  # colorbar 位置 (left, bottom, width, height)
    figsize=(16, 8),               # 图表尺寸
    dendrogram_ratio=0.1,          # 聚类树高度占比
    cbar_kws={'label': 'Expression'}  # colorbar 标签
)

# 🖋️ 设置标题和坐标轴标签（通过 ax_row_dendrogram 和 ax_col_dendrogram 获取）
g.ax_row_dendrogram.set_visible(True)
g.ax_col_dendrogram.set_visible(True)
g.ax_heatmap.set_title('Clustered Heatmap of Combined Top Transcription Factors Across Tissues',
                       fontsize=16, pad=20)
g.ax_heatmap.set_xlabel('Transcription Factor', fontsize=14)
g.ax_heatmap.set_ylabel('Tissue', fontsize=14)

# 🌟 可选：旋转 x 轴标签，避免重叠
plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha='right')

# 📸 展示图表
plt.tight_layout()

plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/TF_heatmap_Tissue-.png')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/TF_heatmap_Tissue-.pdf')

# 构建折线图所需数据：每个 leiden 的 Top1 TF 表达值
# 折线图：展示每个 leiden 的 Top1 TF 在所有 leiden 中的表达趋势
# 指定要展示的 TF 列表
# 假设 adata 已加载，并且你已定义 custom_leiden_order 和 selected_tfs 列表
selected_tfs = ['NME2', 'JUND', 'FOS', 'DLX5']
custom_leiden_order = ['tumor_middle', 'tumor_edge', 'normal_adjacent', 'tumor_metastasis']

# Step 1: 计算每个 leiden 中指定 TF 的平均表达量
mean_exp = []
for leiden in custom_leiden_order:
    cells_in_leiden = adata[adata.obs['leiden'] == leiden]
    mean_per_gene = cells_in_leiden[:, selected_tfs].to_df().mean()
    mean_exp.append(mean_per_gene)

# 构建新的 line_df，行是 leiden，列是每个 TF 的平均表达值
line_df_mean = pd.DataFrame(mean_exp, index=custom_leiden_order, columns=selected_tfs)

# Step 2: 绘图设置
palette = sns.color_palette("husl", n_colors=len(selected_tfs))

plt.figure(figsize=(10, 6))

# 绘制每条折线，使用 TF 名作为 label
for i, tf in enumerate(selected_tfs):
    sns.lineplot(
        x=line_df_mean.index,
        y=line_df_mean[tf],
        marker='o',
        markersize=8,
        linewidth=2,
        color=palette[i],
        label=tf
    )

# 设置坐标轴和标题
plt.xticks(ticks=range(len(custom_leiden_order)), labels=custom_leiden_order, rotation=45, fontsize=12)
plt.xlabel('leiden', fontsize=14)
plt.ylabel('Mean Expression', fontsize=14)  # 改为 Mean
plt.title('Mean Expression Trend of Selected Transcription Factors Across leidens', fontsize=16)

# 美化网格
plt.grid(True, linestyle='--', alpha=0.5)

# 图例移到图外右侧，避免遮挡
plt.legend(title='Transcription Factor', bbox_to_anchor=(1.05, 1), loc='upper left')

# 布局自适应，防止裁剪
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/TF_lineplot_leiden_mean.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/TF_lineplot_leiden_mean.pdf', dpi=300, bbox_inches='tight')
plt.show()
# 美化网格
plt.grid(True, linestyle='--', alpha=0.5)

# 图例移到图外右侧，避免遮挡
plt.legend(title='Transcription Factor', bbox_to_anchor=(1.05, 1), loc='upper left')

# 布局自适应，防止裁剪
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/TF_lineplot_leiden_median.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/TF_lineplot_leiden_median.pdf', dpi=300, bbox_inches='tight')

import scanpy as sc
import pandas as pd
# 添加一列用于分组比较
adata.obs['group'] = adata.obs['leiden'].apply(
    lambda x: 'tumor_edge' if x == 'tumor_edge' else
              ('tumor_middle' if x == 'tumor_middle' else 'others')
)

# 只保留 tumor_edge 和 tumor_middle 进行比较
adata_sub = adata[adata.obs['group'] != 'others'].copy()
# 获取所有 TF 名称（假设你已有一个 TF 列表）
tf_list = transcription_factors  # 示例列表，替换为你实际的 TF 列表

# 限制只分析这些 TF 基因
adata_tf = adata_sub[:, adata_sub.var_names.isin(tf_list)].copy()

# 差异分析
sc.tl.rank_genes_groups(adata_tf, groupby='group', method='wilcoxon', key_added='rank_genes')

# 将结果转换为 DataFrame
result_df = sc.get.rank_genes_groups_df(adata_tf, group='tumor_edge', key='rank_genes')
result_df = result_df[result_df['pvals_adj'] < 0.05]  # 只保留显著差异的 TF
print(result_df.sort_values(by='logfoldchanges', ascending=False))
# 假设 result_df 是 sc.get.rank_genes_groups_df 的结果
sig_up_tfs = result_df[(result_df['pvals_adj'] < 0.05) & (result_df['logfoldchanges'] > 1)]['names'].tolist()
# 提取 tumor_edge 的细胞
tumor_edge_cells = adata_tf[adata_tf.obs['leiden'] == 'tumor_edge']

# 计算每个 TF 的中位表达量
tf_median_in_edge = tumor_edge_cells.to_df().median().loc[sig_up_tfs]
if not tf_median_in_edge.empty:
    highest_exp_tf = tf_median_in_edge.idxmax()
    highest_exp_value = tf_median_in_edge.max()
    print(f"在 tumor_edge 中表达量最高的显著上调 TF 是：{highest_exp_tf}（中位表达值为 {highest_exp_value:.2f}）")
else:
    print("没有满足条件的显著上调 TF。")
