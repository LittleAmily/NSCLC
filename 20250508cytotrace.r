devtools::install_github('zktuong/ktplots', dependencies = TRUE)
.libPaths('~/R/x86_64-pc-linux-gnu-library/4.2')
library(CytoTRACE2)
library(tidyverse)
library(Seurat)
library(paletteer)
library(BiocParallel)
options(bitmapType='cairo')
DimPlot(epi_counts,pt.size = 0.8,group.by = "Detailed_Epithelium_type",label = T)
# 修改Peng_Zhang开头的列名（去除最后两个字符）
peng_cols <- grep("^Peng_Zhang_", colnames(epi_counts), value = TRUE)
new_cols <- substr(peng_cols, 1, nchar(peng_cols) - 2)
colnames(epi_counts)[colnames(epi_counts) %in% peng_cols] <- new_cols
# 读取pathtype.csv文件（注意跳过首行空列名）
pathtype_df <- read.csv("/home/data/sdzl14/NSCLC/zong/epi/pathtype.csv", 
                        row.names = 1, 
                        header = TRUE,
                        check.names = FALSE)
# 假设metadata中有样本名列"sample_name"
# 先初始化Pathtype列为NA
metadata$Pathtype <- NA
# 读取 pathtype.csv 文件
pathtype_file <- read.csv("/home/data/sdzl14/NSCLC/zong/epi/pathtype.csv", header = TRUE, row.names = 1)

# 提取 epi_counts 中以 Peng_Zhang_ 开头的列名
peng_zhang_columns <- grep("^Peng_Zhang_", colnames(epi_counts), value = TRUE)

# 去掉每个列名最后的 "-数字"
processed_column_names <- sub("-\\d+$", "", peng_zhang_columns)

# 创建映射关系: 将处理后的列名与 pathtype 匹配
pathtype_mapping <- sapply(processed_column_names, function(x) {
  # 查找 pathtype.csv 中匹配的行名
  if (x %in% rownames(pathtype_file)) {
    return(pathtype_file[x, "Pathtype"])
  } else {
    return(NA_character_)
  }
})

# 构建 metadata 数据框，保存 pathtype 信息
metadata <- data.frame(
  Sample_ID = processed_column_names,
  Pathtype = pathtype_mapping,
  stringsAsFactors = FALSE
)
# 确保 metadata 的行顺序与 Seurat 对象的列顺序一致（通过 Sample_ID 匹配）
matched_index <- match(colnames(epi_counts), metadata$Sample_ID)

# 将 metadata 中的 Pathtype 信息提取到临时变量
pathtype_from_metadata <- metadata$Pathtype[matched_index]

# 填补 NA 的逻辑（仅替换原数据中 NA 的部分）
epi_counts[['Pathtype']][is.na(epi_counts[['Pathtype']])] <- pathtype_from_metadata[is.na(epi_counts[['Pathtype']])]
saveRDS(epi_counts,file = '~/NSCLC/zong/epi/epi_counts.rds')
lusc <- subset(epi_counts, cells= rownames(epi_counts@meta.data[epi_counts@meta.data$Pathtype=="LUSC",]))
DimPlot(epi_counts,group.by = 'Pathtype')
lusc <- NormalizeData(lusc)
?FindVariableFeatures
lusc <- FindVariableFeatures(lusc,)
lusc <- ScaleData(lusc)
?FindNeighbors
lusc <- FindNeighbors(lusc,reduction = "scanvi_fix_linear")
colnames(epi_counts)
lusc <- FindClusters(lusc,resolution = 0.5)
?RunUMAP

lusc <- RunUMAP(lusc,dims = 1:10)
DimPlot(lusc, reduction = 'umap')
luad <- subset(epi_counts, cells= rownames(epi_counts@meta.data[epi_counts@meta.data$Pathtype=="LUAD",]))

DimPlot(lusc,pt.size = 0.8,group.by = "Detailed_Epithelium_type",label = T)
cytotrace2_res <- cytotrace2(lusc, #seurat对象
                             is_seurat = TRUE, 
                             slot_type = "counts", #counts和data都可以
                             species = 'human')#物种要选择，默认是小鼠
class(cytotrace2_res)
expression_matrix <- as.matrix(epi_counts@assays$RNA@counts)
cytotrace2_result <- cytotrace2(expression_matrix)
cytotrace2_sce <- cytotrace2(epi_counts, 
                             is_seurat = TRUE, 
                             slot_type = "counts", 
                             species = 'human', 
                             seed = 1234)
Idents(cytotrace2_sce)<- 'RNA_snn_res.0.5'
subset <- subset(cytotrace2_sce, cytotrace2_sce$CytoTRACE2_Potency !='Differentiated')
subset = cytotrace2_sce[,cytotrace2_sce@meta.data$CytoTRACE2_Potency %in% c('Unipotent', 'Oligopotent', 'Multipotent')]
DimPlot(subset, group.by = "CytoTRACE2_Potency") 
cluster_0_ <- FindMarkers(cytotrace2_sce,ident.1 = '0',max.cells.per.ident=1000)
comparisons <- combn(unique(cancer@meta.data$RNA_snn_res.1.25), 35, simplify = FALSE)
comparisons
# plotting-一次性生成多个图，然后储存在一个list，用$查看即可
annotation <- data.frame(phenotype = lusc@meta.data$Detailed_Epithelium_type) %>% 
  set_rownames(., colnames(lusc))
plots <- plotData(cytotrace2_result = cytotrace2_res, 
                  annotation = annotation, 
                  is_seurat = TRUE)


plots$CytoTRACE2_UMAP
ggsave("~/NSCLC/zong/fig/MALIGNANT/CytoTRACE2_LUSCUMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Potency_UMAP
ggsave("~/NSCLC/zong/fig/MALIGNANT/CytoTRACE2_Potency_LUSCUMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Relative_UMAP
ggsave("~/NSCLC/zong/fig/MALIGNANT/CytoTRACE2_Relative_LUSCUMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$Phenotype_UMAP
ggsave("~/NSCLC/zong/fig/MALIGNANT/CytoTRACE2_Phenotype_LUSCUMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Boxplot_byPheno
ggsave("~/NSCLC/zong/fig/MALIGNANT/CytoTRACE2_LUSCBoxplot_byPheno.pdf",  width = 9, height = 7, dpi = 300)
luad <- subset(epi_counts, cells= rownames(epi_counts@meta.data[epi_counts@meta.data$Pathtype=="LUAD",]))
DimPlot(epi_counts,group.by = 'Pathtype')
luad <- NormalizeData(luad)
?FindVariableFeatures
luad <- FindVariableFeatures(luad,)
luad <- ScaleData(luad)
?FindNeighbors
luad <- FindNeighbors(luad,reduction = "scanvi_fix_linear")

luad <- FindClusters(luad,resolution = 0.5)
?RunUMAP
luad <- RunPCA(luad)
luad <- RunUMAP(luad,dims = 1:10)
DimPlot(luad, reduction = 'umap')

DimPlot(luad,group.by = 'Detailed_Epithelium_type')
DimPlot(luad,pt.size = 0.8,group.by = "Detailed_Epithelium_type",label = T)

subset <- subset(luad,downsample = 500)
DimPlot(subset,pt.size = 0.8,group.by = "Detailed_Epithelium_type",label = T)
cytotrace2_res <- cytotrace2(luad, #seurat对象
                             is_seurat = TRUE, 
                             slot_type = "counts", #counts和data都可以
                             species = 'human')#物种要选择，默认是小鼠

# plotting-一次性生成多个图，然后储存在一个list，用$查看即可
annotation <- data.frame(phenotype = luad@meta.data$Origin) %>% 
  set_rownames(., colnames(luad))
plots <- plotData(cytotrace2_result = cytotrace2_res, 
                  annotation = annotation, 
                  is_seurat = TRUE)

luad@meta.data[["CytoTRACE2_Score"]] <- cytotrace2_res@meta.data[["CytoTRACE2_Score"]]
luad@meta.data[["CytoTRACE2_Potency"]] <- cytotrace2_res@meta.data[["CytoTRACE2_Potency"]]
luad@meta.data[["CytoTRACE2_Relative"]] <- cytotrace2_res@meta.data[["CytoTRACE2_Relative"]]

# 使用 FeaturePlot 绘制 CytoTRACE2_Relative 特征，并调整颜色渐变
p <- FeaturePlot(
  object = luad,
  features = "CytoTRACE2_Score",
  cols = viridis::viridis(3),  # 使用 viridis 调色板
  pt.size = 0.5,                  # 设置点的大小
  raster = T                  # 是否使用栅格化绘图
)

# 显示图形
print(p)
plots$CytoTRACE2_UMAP
ggsave("~/NSCLC/zong/fig/MALIGNANT/CytoTRACE2_luadUMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Potency_UMAP
ggsave("~/NSCLC/zong/fig/MALIGNANT/CytoTRACE2_Potency_luadUMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Relative_UMAP
ggsave("~/NSCLC/zong/fig/MALIGNANT/CytoTRACE2_Relative_luadUMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$Phenotype_UMAP
ggsave("~/NSCLC/zong/fig/MALIGNANT/CytoTRACE2_Phenotype_luadUMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Boxplot_byPheno
ggsave("~/NSCLC/zong/fig/MALIGNANT/CytoTRACE2_luadBoxplot_byPheno.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2
plots$CytoTRACE2_Potency
plots$CytoTRACE2_Relative
plots$CytoTRACE2_Relative_style 
plots$Phenotype
plots$CytoTRACE2_Boxplot_byPheno
library(harmony)
sce_harmony <-  cancer %>% RunHarmony("orig.ident")
DimPlot(sce_harmony, group.by = "PatientID", reduction = 'harmony')
cytotrace2_sce1 <- cytotrace2(sce_harmony, 
                              is_seurat = TRUE, 
                              slot_type = "counts", 
                              species = 'human', 
                              seed = 1234)
plots1 <- CytoTRACE2_plotData( seurat = sce_harmony, # single cell object
                               reduction = "harmony", #reduction改为harmony
                               cytotrace2_result = cytotrace2_sce1,  # cytotrace2 object
                               group = 'RNA_snn_res.0.5',  # 分组
                               add_jitter =F, #不加抖动点
                               is_seurat = TRUE)
plots1$CytoTRACE2
plots1$CytoTRACE2_Potency
plots1$CytoTRACE2_Relative
plots1$CytoTRACE2_Relative_style 
plots1$Phenotype
plots1$CytoTRACE2_Boxplot_byPheno
save(cancer_14968_harmony,cluster_0,cluster_1,cluster_11,cluster_2,cluster_3,cluster_4,cluster_7,cluster_9,file = '~/NSCLC/scRNA_data_count/cancer分群marker.rda')
#无脑运行即可，无需解释
source(file = '~/NSCLC/Peilin_Wang_2025/code/Vector.R')
options(bitmapType='cairo')
.libPaths('~/R/x86_64-pc-linux-gnu-library/4.2')
?install_local
devtools::install_local('~/NSCLC/scRNA_data_count/rigraph-main/rigraph-1.6.0')
options(bitmapType='cairo')
VEC = epi_counts@reductions$umap@cell.embeddings #提取umap坐标
rownames(VEC) = colnames(epi_counts)
PCA = epi_counts@reductions$pca@cell.embeddings
# Remove quantile-based colinearity among PCs (new feature in VECTOR 0.0.3):   
PCA=vector.rankPCA(PCA)
# Define pixel
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
# Build network
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
# Calculate Quantile Polarization (QP) score
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
# Get pixel's QP score
OUT=vector.gridValue(OUT,SHOW=TRUE)
# Find starting point
V(OUT$GRAPH)$name
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
# Infer vector
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
