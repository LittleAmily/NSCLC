import scanpy as sc
import pandas as pd
import gseapy as gp
from gseapy import prerank, barplot, dotplot

# Step 1: 加载 AnnData 文件
adata = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/malignant_integrated.h5ad')
adata = adata.copy()
# Step 2: 获取表达矩阵（这里以平均表达为例）
# Step 2: 过滤只保留指定 Origin
print(adata.obs['Origin'].cat.categories)
custom_Origin_order = ['adrenal', 'brain', 'liver', 'lung', 'lymph_node', 'pleura/effusion']
adata = adata[adata.obs['Origin'].isin(custom_Origin_order)].copy()

# Step 3: 差异分析（使用 Wilcoxon test）
sc.tl.rank_genes_groups(adata, groupby='Origin', method='wilcoxon', key_added='rank_genes_Origin')

# Step 4: 提取每个 Origin 的 log2FoldChange 并构建 rank 列表
rank_dict = {}

for Origin in adata.obs['Origin'].cat.categories:
    # 获取该 Origin 和其他 Origin 的表达值
    df = sc.get.rank_genes_groups_df(adata, group=Origin, key='rank_genes_Origin')
    df.set_index('names', inplace=True)
    rank_dict[Origin] = df['scores']
    # 查看某一个 Origin 的 DEG 结果
df_deg = sc.get.rank_genes_groups_df(adata, group='0', key='rank_genes')
print(df_deg.head())
# scores 即 log2FoldChange# Step 3: 解析 GMT 文件
gmt_path = '/home/data/sdzl14/NSCLC/zong/h.all.v2024.1.Hs.symbols.gmt.txt'

hallmark_gene_sets = {}

with open(gmt_path, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        pathway_name = parts[0]
        genes = parts[2:]  # 第三列开始是基因名
        hallmark_gene_sets[pathway_name] = genes

# Step 4: 执行 prerank GSEA 分析
results = []

for Origin in custom_Origin_order:
    print(f"\n--- Processing Origin: {Origin} ---")
    
    # 获取当前 Origin 的排序基因列表
    rank_list = rank_dict[Origin].reset_index()
    rank_list.columns = ['gene', 'log2FC']
    rank_list = rank_list.drop_duplicates(subset='gene').set_index('gene')

    print(f"{Origin}: Total genes available = {len(rank_list)}")

    # 确保 rank_list 非空
    if rank_list.empty:
        print(f"[ERROR] No valid genes found for {Origin}. Skipping...")
        continue

    # 执行 prerank GSEA（添加 try-except 避免崩溃）
    try:
        pre_res = gp.prerank(
            rnk=rank_list,
            gene_sets=hallmark_gene_sets,
            outdir=None,
            seed=6,
            verbose=True
        )

        # 检查 pre_res 是否为空
        if pre_res is None:
            print(f"[ERROR] prerank returned None for {Origin}")
            continue

        # 提取当前 Origin 的富集结果
        Origin_results = []
        for term in pre_res.results:
            res = pre_res.results[term]
            Origin_results.append([
                Origin,
                term,
                res['fdr'],
                res['es'],
                res['nes']
            ])

        # 转换为 DataFrame 并排序
        Origin_df = pd.DataFrame(Origin_results, columns=['Origin', 'Term', 'fdr', 'es', 'nes'])
        Origin_df = Origin_df.sort_values('fdr')
        results.append(Origin_df)

    except Exception as e:
        print(f"[EXCEPTION] Error in prerank for {Origin}: {str(e)}")
        continue

# Step 6: 合并所有 Origin 的结果
if results:
    final_result = pd.concat(results, axis=0).reset_index(drop=True)
    print("\nFinal Result:")
    print(final_result.head())
else:
    print("[WARNING] No results were generated.")
from gseapy.plot import dotplot
final_result.to_csv('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/gsea_resultsOrigin.csv', index=False)
# Step 1: 转换格式为 dotplot 支持的格式
df = final_result.pivot(index='Term', columns='Origin', values='nes')

# Step 2: 使用 dotplot 可视化
dotplot(df,
        title='Hallmark Pathway Enrichment (Normalized Enrichment Score)',
        cmap='vlag',
        figsize=(8, 10),
        cutoff=0.05,
        ofname='/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/gsea_dotplot.png')
heatmap_data = final_result.pivot(index='Term', columns='Origin', values='nes')
import seaborn as sns
import matplotlib.pyplot as plt
# 假设 final_result 是一个 pandas DataFrame
# 提取所有 Origin 类型
Origins = final_result['Origin'].unique()

# 找出在所有 Origin 中都存在的 Term
common_terms = set(final_result[final_result['Origin'] == Origins[0]]['Term'])

for Origin in Origins[1:]:
    common_terms &= set(final_result[final_result['Origin'] == Origin]['Term'])

# 过滤数据以仅包括 common_terms
filtered_data = final_result[final_result['Term'].isin(common_terms)]

# 创建一个数据透视表用于热图
heatmap_data = filtered_data.groupby(['Term', 'Origin'])['nes'].mean().unstack()
# 使用 clustermap 进行聚类热图绘制，同时更改配色方案为 "coolwarm"
g = sns.clustermap(
    heatmap_data,
    cmap="coolwarm",            # 更好看的配色
    annot=True,                 # 显示数值
    annot_kws={"size": 8},     # 注释字体大小
    linewidths=.5,              # 单元格边框宽度
    figsize=(10, 15),           # 图片大小
    cbar_kws={"shrink": .5},    # 调整颜色条尺寸
    dendrogram_ratio=0.1,        # 控制树状图比例
    yticklabels=True,  
)

# 设置标题和标签
g.ax_row_dendrogram.set_visible(False)  # 隐藏行树状图（可选）
g.ax_col_dendrogram.set_visible(False)  # 隐藏列树状图（可选）

# 调整布局并保存
plt.suptitle('Clustered Heatmap of NES Values for Common Terms', y=1.02)
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/gsea_clustered_heatmapOrigin.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/gsea_clustered_heatmapOrigin.pdf', dpi=300, bbox_inches='tight')

# 绘制热图
plt.figure(figsize=(15, 8))
sns.heatmap(heatmap_data, annot=True, cmap="viridis")
plt.title('Heatmap of NES Values for Common Terms')
plt.xlabel('Terms')
plt.ylabel('Origin')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/gsea_heatmap.png')
# 绘制点图
plt.figure(figsize=(12, 6))
sns.scatterplot(data=filtered_data, x='Term', y='Origin', size='nes', hue='fdr', sizes=(20, 200), legend=False)
plt.title('Dotplot of NES and FDR for Common Terms')
plt.xlabel('Terms')
plt.ylabel('Origin')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/gsea_heatmap.png')
# 感兴趣的通路列表
selected_terms = [
    'HALLMARK_E2F_TARGETS',
    'HALLMARK_MYC_TARGETS_V1',
    'HALLMARK_G2M_CHECKPOINT',
    'HALLMARK_HEDGEHOG_SIGNALING',
    'HALLMARK_MITOTIC_SPINDLE',
    'HALLMARK_UV_RESPONSE_DN',
    'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
    'HALLMARK_COAGULATION',
    'HALLMARK_IL6_JAK_STAT3_SIGNALING',
    'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
    'HALLMARK_NOTCH_SIGNALING',
    'HALLMARK_GLYCOLYSIS'
]

# 筛选 heatmap_data 中这些 Term 的数据
subset_heatmap_data = heatmap_data.loc[selected_terms]

# 绘制 clustermap（可选 dendrogram）
g = sns.clustermap(
    subset_heatmap_data,
    cmap="coolwarm",
    annot=True,
    annot_kws={"size": 8},
    linewidths=.5,
    figsize=(10, 8),           # 根据 Term 数量调整高度
    cbar_kws={"shrink": .5},
    dendrogram_ratio=0.1,
    yticklabels=True,
    row_cluster=False,         # 不聚类行（如果不需要）
    col_cluster=True           # 列聚类保持开启
)

# 设置标题和隐藏树状图（可选）
g.ax_row_dendrogram.set_visible(False)
g.ax_col_dendrogram.set_visible(False)

# 保存图片
plt.suptitle('Selected Hallmark Pathways', y=1.02)
plt.tight_layout()
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/gsea_selected_pathwaysOrigin.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/gsea_selected_pathwaysOrigin.pdf', dpi=300, bbox_inches='tight')
# 筛选 final_result 中这些 Term 的数据
dotplot_data = final_result[final_result['Term'].isin(selected_terms)]

# 绘制点图
plt.figure(figsize=(12, 8))
scatter = sns.scatterplot(
    data=dotplot_data,
    x='Term',
    y='Origin',
    size='nes',
    hue='fdr',
    sizes=(50, 600),               # 调整点的大小范围，避免过大
    palette="RdYlBu_r",            # 更美观的颜色映射（红-黄-蓝反转）
    alpha=0.85,
    edgecolor='black',
    linewidth=0.5
)

# 设置图形细节
plt.title('Dotplot of NES and FDR for Selected Hallmark Pathways', fontsize=16)
plt.xlabel('Pathway (Term)', fontsize=14)
plt.ylabel('Origin Type', fontsize=14)
plt.xticks(rotation=45, ha='right', fontsize=12)
plt.yticks(fontsize=12)

# 调整布局，防止被裁剪
plt.tight_layout()

# 添加图例（分离颜色和大小）
handles, labels = scatter.get_legend_handles_labels()
# 找出 color 和 size 对应的图例项
num_fdr_labels = int(len(labels) / 2)

legend1 = plt.legend(
    handles[:num_fdr_labels], 
    labels[:num_fdr_labels],
    title="FDR", 
    loc="upper left", 
    bbox_to_anchor=(1, 1),
    frameon=False,
    fontsize=12,
    title_fontsize=13
)

legend2 = plt.legend(
    handles[num_fdr_labels:], 
    labels[num_fdr_labels:],
    title="NES", 
    loc="lower left", 
    bbox_to_anchor=(1, 0),
    frameon=False,
    fontsize=12,
    title_fontsize=13
)

plt.gca().add_artist(legend1)

# 保存图像
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/dotplot_selected_pathwaysorigin.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/dotplot_selected_pathwaysorigin.pdf', dpi=300, bbox_inches='tight')

plt.show()
gp.gseaplot(rank_metric = pre_res.ranking,term = 'HALLMARK_WNT_BETA_CATENIN_SIGNALING', **pre_res.results['HALLMARK_WNT_BETA_CATENIN_SIGNALING'],color="black")
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/MALIGNANT/gseaplot.png', dpi=300, bbox_inches='tight')
ad = sc.read_h5ad('/home/data/sdzl14/NSCLC/zong/stromal.h5ad')
ad = ad.copy()
ad.obs['_prediction'].cat.categories
sc.pl.umap(ad, color=["stromal_celltype"])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/stromal/stromal_celltype.png',dpi = 1200, bbox_inches='tight')
# 假设你已经有了 CAF 对象，并且已经定义了 CAF_type

df = ad.obs[['Tissue', 'Epithelium_type']].copy()
ad.obs['Tissue']
# 固定 Tissue 的顺序
desired_order =['tumor_middle', 'tumor_edge', 'normal_adjacent', 'tumor_metastasis']
df['Tissue'] = pd.Categorical(df['Tissue'], categories=desired_order, ordered=True)

# 构建交叉表并计算比例
count_df = pd.crosstab(index=df['Tissue'], columns=df['Epithelium_type'])
prop_df = count_df.div(count_df.sum(axis=1), axis=0)
#创建新的 DataFrame 来存储处理后的比例数据
prop_df_with_other = pd.DataFrame(index=prop_df.index, columns=prop_df.columns.tolist() + ['Other'])

# 遍历每一行（每个 Timepoint "Other"
threshold = 0.10
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
plt.title('Proportion of B Subtypes by Tissue (Specified Order)', fontsize=14)
plt.xlabel('Tissue')
plt.ylabel('Proportion')
plt.xticks(rotation=45)
plt.legend(title='B Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()



# 保存图像
results_file = '/home/data/sdzl14/NSCLC/zong/fig/EPI/'
plt.savefig(results_file + 'EPi_Tissue_immune_type_stacked_barplot_ordered.png', dpi=300, bbox_inches='tight')
plt.savefig(results_file + 'EPi_Tissue_immune_type_stacked_barplot_ordered.pdf', dpi=300, bbox_inches='tight')






plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14 
sc.pl.umap(ad,color=["Origin"])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/EPI/Origin-.png',dpi = 1200,bbox_inches = 'tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/QC/Origin-.pdf',dpi = 300,bbox_inches = 'tight')
adata = ad[ad.obs['Tissue'] == 'normal_distant']
sc.pl.umap(adata,color=["n_genes_by_counts"])
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/QC/UMAPnormal_distant.png',dpi = 1200,bbox_inches = 'tight')
plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/QC/UMAPnormal_distant.pdf',dpi = 300,bbox_inches = 'tight')
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

# 随机抽样10000个细胞
if ad.n_obs > 100000:
    # 使用随机抽样
    sc.pp.subsample(ad, n_obs=100000, random_state=123)
else:
    print(f"数据集只有 {ad.n_obs} 个细胞，使用全部细胞")

# 执行差异分析
sc.tl.rank_genes_groups(ad, 'Celltype', method='wilcoxon', key_added='rank_genes_celltype')
ad.obs['Celltype'] = ad.obs['Celltype'].astype('category')
ad.obs['Celltype'].cat.categories
# 提取每个细胞类型的top1基因
top_genes = []
for celltype in ad.obs['Celltype'].cat.categories:
    df = sc.get.rank_genes_groups_df(ad, group=celltype, key='rank_genes_celltype')
    
    # 跳过空结果
    if df.empty:
        print(f"警告: {celltype} 组没有差异基因")
        continue
    
    # 按logfoldchange排序并取top1
    top_gene = df.sort_values('scores', ascending=False).iloc[0]['names']
    top_genes.append(top_gene)

# 去重并过滤空值
top_genes = list(np.unique([g for g in top_genes if g]))

# 检查基因是否在数据中
missing_genes = [g for g in top_genes if g not in ad.var_names]
if missing_genes:
    print(f"警告: 以下基因不在数据中: {', '.join(missing_genes)}")
    top_genes = [g for g in top_genes if g in ad.var_names]


# 绘制热图（添加错误处理）
try:
    sc.pl.heatmap(
        ad,
        var_names=top_genes,
        groupby='Celltype',
        show_gene_labels=True,
        figsize=(12, 8),
        cmap='viridis',
        dendrogram=True,
        swap_axes=False,
        use_raw=False,
        title='Top 1 DEG per Celltype',
        vmin=-2,
        vmax=2
    )
    plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/top1_genes_heatmap.png', dpi=300, bbox_inches='tight')
    
except Exception as e:
    print(f"绘制热图时出错: {str(e)}")
    # 备选方案：使用矩阵热图
    plt.figure(figsize=(12, 8))
    sc.pl.matrixplot(
        ad,
        var_names=top_genes,
        groupby='Celltype',
        cmap='viridis',
        title='Top 1 DEG per Celltype'
    )
    plt.savefig('/home/data/sdzl14/NSCLC/zong/fig/top1_genes_matrixplot.png', dpi=300, bbox_inches='tight')
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list

# 1. 创建表达矩阵
# 获取每个细胞类型的平均表达值（使用log1p标准化后的数据）
expr_matrix = pd.DataFrame(
    ad.X.toarray() if hasattr(ad.X, 'toarray') else ad.X,
    index=ad.obs['Celltype'],
    columns=ad.var_names
)[top_genes].groupby(level=0).mean()

# 2. 确保矩阵行和列的顺序对齐（使对角线高表达）
# 创建基因到细胞类型的映射
gene_to_celltype = {}
for celltype in ad.obs['Celltype'].cat.categories:
    # 获取该细胞类型的top基因
    df = sc.get.rank_genes_groups_df(ad, group=celltype, key='rank_genes_celltype')
    if not df.empty:
        top_gene = df.sort_values('logfoldchanges', ascending=False).iloc[0]['names']
        gene_to_celltype[top_gene] = celltype

# 3. 重新排序矩阵，使每个基因与其所属细胞类型对齐
# 获取唯一的细胞类型列表
celltypes = ad.obs['Celltype'].cat.categories.tolist()

# 为每个细胞类型选择其top基因
diagonal_genes = []
for ct in celltypes:
    # 找到属于该细胞类型的top基因
    genes_for_ct = [gene for gene, ct_map in gene_to_celltype.items() if ct_map == ct]
    if genes_for_ct:
        diagonal_genes.append(genes_for_ct[0])  # 取第一个（即top1）

# 过滤矩阵，只包含对角基因
diagonal_matrix = expr_matrix.loc[celltypes, diagonal_genes]

# 4. 绘制优化的聚类热图
plt.figure(figsize=(14, 12))

# 创建自定义行和列链接
row_linkage = linkage(diagonal_matrix, method='average', metric='euclidean')
col_linkage = linkage(diagonal_matrix.T, method='average', metric='euclidean')

# 获取聚类顺序
row_order = leaves_list(row_linkage)
col_order = leaves_list(col_linkage)

# 重新排序矩阵
ordered_matrix = diagonal_matrix.iloc[row_order, col_order]

# 绘制热图
sns.heatmap(
    ordered_matrix,
    cmap='viridis',
    annot=True,
    fmt=".1f",
    annot_kws={"size": 8},
    linewidths=0.5,
    linecolor='lightgray',
    cbar_kws={"label": "Expression Level (z-score)"}
)

# 5. 添加对角线标记
for i, (celltype, gene) in enumerate(zip(ordered_matrix.index, ordered_matrix.columns)):
    # 标记对角线单元格
    plt.text(
        i + 0.5, i + 0.5, 
        f"{gene}\n({celltype})", 
        ha='center', va='center',
        fontsize=9, color='white', 
        bbox=dict(facecolor='black', alpha=0.5, boxstyle='round,pad=0.3')
    )

# 6. 设置标签和标题
plt.title('Top 1 DEG per Celltype (Diagonal Clustermap)', fontsize=16, pad=20)
plt.xlabel('Genes', fontsize=14)
plt.ylabel('Cell Types', fontsize=14)
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.yticks(fontsize=10)

# 7. 保存高质量图像
plt.tight_layout()
plt.savefig(
    '/home/data/sdzl14/NSCLC/zong/fig/top1_genes_diagonal_clustermap.png', 
    dpi=300, 
    bbox_inches='tight',
    transparent=True
)
plt.show()
