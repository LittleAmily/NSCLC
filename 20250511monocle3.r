setwd('~/NSCLC/')
devtools::install_github('cole-trapnell-lab/monocle3')
.libPaths('~/R/x86_64-pc-linux-gnu-library/4.2')
library(monocle3)
library(tidyverse)
library(reticulate)
library(CytoTRACE)
require(Seurat)
require(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
library(tidydr)
library(CytoTRACE2)
library(SeuratData)
library(remotes)
library(ggsci)
library(matrixStats)
epi_counts <- readRDS("~/NSCLC/zong/epi/epi_counts.rds")
lusc <- subset(epi_counts, cells= rownames(epi_counts@meta.data[epi_counts@meta.data$Pathtype=="LUSC",]))

lusc <- subset(epi_counts,downsample = 30000)
expression_matrix = lusc@assays$RNA@data
cell_metadata = data.frame(lusc@meta.data)
gene_annotation = data.frame(expression_matrix[,1])
gene_annotation[,1] = row.names(gene_annotation)
##构建Monocle3 cds对象
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#data数据选择norm_method = c("none")
cds <- preprocess_cds(cds, num_dim = 50,norm_method = c("none"))
plot_pc_variance_explained(cds)
#去除批次效应,多个样本时可以通过此方式去除批次
cds <- align_cds(cds, alignment_group = "Sample")
?reduce_dimension
## 降维，默认是"Umap"方式
cds <- reduce_dimension(cds,cores=20,umap.min_dist = 0.3,umap.n_neighbors = 30)
plot_cells(cds)
## 聚类分群，分辨率调小，是为了让细胞是一群可以更好展示
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
?cluster_cells
##这里我们假定以B细胞为起点
myselect <- function(cds,select.classify,my_select){
  cell_ids <- which(colData(cds)[,select.classify] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes}

cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'Detailed_Epithelium_type',my_select = "Basal cell"))
##拟时序值越高表示细胞分化程度越高，这里仅为演示，并非真实分化情况
plot_cells(cds, color_cells_by = "pseudotime",
           show_trajectory_graph=F)
ggsave(plot_cells(cds, color_cells_by = "pseudotime",
                  show_trajectory_graph=F),file = '~/NSCLC/zong/fig/EPI/monocle3_luad-.pdf',height = 6,width = 6)
ggsave(plot_cells(cds,
                                                 color_cells_by = "Timepoint",
                                                 label_cell_groups=FALSE,
                                                 label_leaves=FALSE,
                                                 label_branch_points=FALSE,
                                                 graph_label_size=0.1,
                                                 show_trajectory_graph=F),file = '~/NSCLC/zong/fig/EPI/monocle3_luadTimepoint.pdf',height = 6,width = 6)                 
plot_cells(cds, color_cells_by="cluster_cells")

colnames(gene_annotation)=c("gene_short_name")
DimPlot(epi_counts,label = T,
        cols = c(
          ,pal_d3('category20b',alpha = 0.6)(15))) 
expression_matrix = epi_counts@assays$RNA@data
cell_metadata = data.frame(epi_counts@meta.data)
gene_annotation = data.frame(expression_matrix[,1])
gene_annotation[,1] = row.names(gene_annotation)
colnames(gene_annotation)=c("gene_short_name")
##构建Monocle3 cds对象
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#data数据选择norm_method = c("none")
cds <- preprocess_cds(cds, num_dim = 50,norm_method = c("none"))
plot_pc_variance_explained(cds)
#去除批次效应,多个样本时可以通过此方式去除批次
cds <- align_cds(cds, alignment_group = "Sample")

## 降维，默认是"Umap"方式
cds <- reduce_dimension(cds,cores=20)
plot_cells(cds)
## 聚类分群，分辨率调小，是为了让细胞是一群可以更好展示
cds <- cluster_cells(cds)
?cluster_cells
plot_cells(cds, color_cells_by="cluster_cells")
?learn_graph
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
## 拟时序
cds <- learn_graph(cds)
save(cancer,cds,file = 'monocle3.rda')
##选择特定细胞作为起点，代码比交互式页面方便许多，且不会出现选中自己不想要的细胞
##这里我们假定以B细胞为起点
myselect <- function(cds,select.classify,my_select){
  cell_ids <- which(colData(cds)[,select.classify] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes}

cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'Detailed_Epithelium_type',my_select = "Basal cell"))
##使用Seurat的UMAP信息，这样可以与Seurat对象的细胞分布保持一致
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(epi_counts, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

##不同细胞类型拟时序数值
##拟时序值越高表示细胞分化程度越高，这里仅为演示，并非真实分化情况
plot_cells(cds, color_cells_by = "pseudotime",
           show_trajectory_graph=F) + plot_cells(cds,
                                                 color_cells_by = "Detailed_Epithelium_type",
                                                 label_cell_groups=FALSE,
                                                 label_leaves=FALSE,
                                                 label_branch_points=FALSE,
                                                 graph_label_size=0.1,
                                                 show_trajectory_graph=F)
#差异基因展示
Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=6)
Track_genes$gene_short_name <-rownames(Track_genes)
#按莫兰指数选择TOP基因
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>%pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by="RNA_snn_res.1.25",min_expr=0.5, ncol= 2,cell_size=1.5) + 
  scale_color_manual(values = c(pal_jco("default",alpha = 0.6)(10),pal_d3("category20b",alpha = 0.6)(20),pal_d3('category20b',alpha = 0.6)(2))) 
Track_genes_sig ='GPNMB'
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=TRUE,  label_leaves=FALSE)
ciliated_genes <- c("GPNMB",
                    "SET",
                    "IGLC2",
                    "NME2",
                    "ALDOA",
                    "CLU",
                    "S100A2",
                    "S100A8",
                    "S100A9",
                    "S100A10",
                    "YBX1",
                    "PTPRC","PPDPF","FTL")
cds_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime")
fit_coefs <- coefficient_table(gene_fits)
fit_coefs <- coefficient_table(gene_fits, component = "status")
pseudotime_matrix <- cds@principal_graph_aux$UMAP$pseudotime
# 提取 cds 中的细胞ID（Monocle3 顺序）
cds_cell_ids <- colnames(cds)

# 提取 Seurat 对象中的细胞ID（lusc 子集）
lusc_cell_ids <- colnames(lusc)

# 验证细胞ID是否完全匹配（顺序可能不同）
identical(sort(cds_cell_ids), sort(lusc_cell_ids))  # 应返回 TRUE

# 若上述返回 TRUE，直接匹配
lusc$pseudotime <- pseudotime_matrix[match(lusc_cell_ids, cds_cell_ids)]
# 查看元数据
head(lusc@meta.data[, c("Tissue", "pseudotime")])  # 替换 Cell_ID_column 为实际列名

# 提取拟时序向量（确保已正确添加到 lusc 对象）
pseudotime <- lusc$pseudotime

# 提取目标基因表达矩阵（log标准化后的数据）
gene_list <- c("GPNMB","SET","IGLC2","NME2","ALDOA","CLU","S100A2",
               "S100A8","S100A9","S100A10","YBX1","PTPRC","PPDPF","FTL")
               # 方法一：直接从scale.data提取（需确认已执行过ScaleData）
library(Seurat)
lusc <- NormalizeData(lusc)
lusc <- ScaleData(lusc)
expression_scaled <- lusc@assays[["RNA"]]@scale.data
# 提取目标基因数据
expression_data <- expression_scaled[intersect(gene_list, rownames(expression_scaled)), ]
rownames(expression_data)
expression_data <- t(expression_data)
# 创建分析数据框
analysis_df <- data.frame(
  Pseudotime = pseudotime,
  expression_data
) %>%
  tidyr::pivot_longer(
    cols = -Pseudotime,
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  filter(!is.na(Pseudotime)) 

# 步骤 1: 数据清洗
# 检查Expression列的潜在嵌套结构
if(any(map_lgl(analysis_df$Expression, ~!is.atomic(.x)))) {
  message("检测到非原子类型数据，需展开嵌套结构")
  analysis_df <- analysis_df %>% 
    tidyr::unnest(Expression) %>%  # 明确指定展开Expression列
    tidyr::unnest(Pseudotime)      # 确保所有列都是展开状态
}
analysis_clean <- analysis_df %>%
  mutate(
    Expression = map_dbl(Expression, ~{
      val <- suppressWarnings(as.numeric(.x))
      ifelse(is.na(val), NA_real_, val)  # 保留NA用于后续过滤
    })
  ) 

# 最终安全版本
analysis_clean <- analysis_df %>%
  filter(!is.infinite(Pseudotime))
head(analysis_clean)
# 步骤 2: 安全拟合

fit_models <- analysis_clean %>%
  group_by(Gene) %>%
  do({
    df <- .
    tryCatch({
      mod <- gam(
        Expression ~ s(Pseudotime, k = 3, bs = "tp"),  # 显式指定基函数
        data = df,
        method = "REML"
      )
      tidy(mod) %>%
        filter(term == "s(Pseudotime)") %>%
        mutate(p.value = approx.p.value)
    }, error = function(e) {
      data.frame(term = NA, p.value = NA)
    })
  }) %>%
  filter(!is.na(term))# 步骤 3: 安全拟合模型
head(fit_models)
# 过滤 Inf 值（假设所有 Inf 为无效值）
analysis_clean <- analysis_df %>%
  filter(
    !is.na(Pseudotime),
    !is.infinite(Pseudotime),
    Pseudotime > 0  # 可选：去除零或负值（根据实际分布调整）
  )
# 绘制p值的柱状图
# 定义一个函数来计算 R² 和 p 值
calculate_stats <- function(df) {
  mod <- gam(Expression ~ s(Pseudotime, k = 3, bs = "tp"), data = df, method = "REML")
  summary_mod <- summary(mod)
  
  # 提取 R²
  r_squared <- summary_mod$r.sq
  
  # 提取 p 值
  p_value <- summary_mod$s.table["s(Pseudotime)", "p-value"]
  
  return(data.frame(r_squared = r_squared, p_value = p_value))
}

# 计算每个基因的统计信息
gene_stats <- analysis_clean %>%
  group_by(Gene) %>%
  do(calculate_stats(.))

# 合并统计信息到原始数据框
analysis_clean_with_stats <- left_join(analysis_clean, gene_stats, by = "Gene")
head(analysis_clean_with_stats)
# 绘制每个基因与 Pseudotime 的相关性图片，并添加 R² 和 p 值
plot_gene_correlation <- function(data, gene_name) {
  sub_data <- data %>% filter(Gene == gene_name)
  
  # 提取 r_squared 和 p_value
  r2 <- unique(sub_data$r_squared)
  p_val <- unique(sub_data$p_value)
  
  # 格式化标签
  label <- paste0("R² = ", round(r2, 3), "\np = ", format.pval(p_val, digits = 2))
  
  ggplot(sub_data, aes(x = Pseudotime, y = Expression)) +
    geom_point(aes(color = Expression), alpha = 0.7, size = 2) + # 使用表达值作为颜色
    scale_color_gradient(low = "blue", high = "red") + # 设置颜色渐变
    geom_smooth(method = "lm", color = "black", se = FALSE, size = 0.5) + # 黑色回归线
    labs(title = gene_name,
         x = "Pseudotime",
         y = "Expression",
         color = "Expression Level") +
    annotate("text", x = Inf, y = -Inf, label = label, hjust = 1, vjust = -0.5,
             parse = FALSE, family = "mono") +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey90"), # 调整网格线颜色
          panel.grid.minor = element_blank(), # 移除次要网格线
          text = element_text(size = 10), # 统一文本大小
          plot.title = element_text(hjust = 0.5)) # 居中标题
}
genes_to_plot <- c("GPNMB","SET","IGLC2","NME2","ALDOA","CLU","S100A2",
                   "S100A8","S100A9","S100A10","YBX1","PTPRC","PPDPF","FTL")

# 使用 lapply 生成每个基因的图
plots <- lapply(genes_to_plot, function(gene) plot_gene_correlation(analysis_clean_with_stats, gene))

# 拼接成一张大图（每行 4 张图）
final_plot <- wrap_plots(plots, ncol = 4)
ggsave("~/NSCLC/zong/fig/EPI/gene_pseudotime_correlation.pdf", plot = final_plot, width = 20, height = 15)
x# 显示图像
print(final_plot)
# 使用 png 设备保存图像
png("~/NSCLC/zong/fig/EPI/Gene_vs_Pseudotime_Correlation.png", width = 10.1, height = 10.1, units = "in", res = 300)
print(p)
dev.off()
options(bitmapType='cairo')
