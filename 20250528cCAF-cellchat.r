setwd('~/NSCLC/')
options(bitmapType='cairo')
.libPaths('~/R/x86_64-pc-linux-gnu-library/4.2')
.libPaths()
devtools::install_github('cole-trapnell-lab/monocle3')
# 卸载当前 Matrix 版本
remove.packages("Matrix")
install.packages("Matrix")
install.packages('/home/sdzl/jinworks-CellChat-v2.1.2-5-g346fb61.tar.gz', repos = NULL, type = "source")
# 安装旧版 Matrix
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-0.tar.gz", 
                 repos = NULL, type = "source")
library(Seurat)
library(CellChat)
library(tidyverse)
library(patchwork)
devtools::install_github("jinworks/CellChat")
library(scPRIT)
library(ggsci)
devtools::install_local("/Users/zhangjinyang/Desktop/CellChat-main.zip")#后面这个就是包的路径
devtools::install_local("~/NSCLC/CellChat-main-2.zip")
library(svglite)
tumor <- readRDS('/home/data/sdzl14/NSCLC/zong/malignant_integrated.rds')
tumor@meta.data$celltype_coarse <- 'Tumor'
tumor@meta.data$immune_celltype <- 'Tumor'
CAF <- readRDS('/home/data/sdzl14/NSCLC/zong/CAF.rds')
CAF@meta.data$celltype_coarse <- CAF@meta.data$CAF_type 
dim(immune)
# 合并Seurat对象
merged_object <- merge(x = tumor, y = CAF,  add.cell.ids = c("Tumor", "CAF"))
Idents(merged_object) <- 'celltype_coarse'
levels(Idents(merged_object) )
NormalizeData(merged_object)
ScaleData(merged_object)
Idents(merged_object) <- 'Tissue'
subset <- subset(merged_object,ident = 'tumor_middle')
data.input <- GetAssayData(subset, slot = 'data') # normalized data matrix
meta <- subset@meta.data[,c("Sample","celltype_coarse","Origin")]
colnames(meta) <-  c("Sample","labels","Origin")
identical(rownames(meta),colnames(data.input))
# 构建cellchat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
levels(cellchat@idents)
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
cellchat@DB <- CellChatDB.use
# Only uses the Secreted Signaling from CellChatDB v1
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole databa
future::plan("multisession", workers = 10) # do parallel
cellchat@data.signaling[cellchat@data.signaling <= 0] <- 0.01  # 加一个 pseudo-count 避免 log(0)
options(future.globals.maxSize = 8 * 1024^3)  # 设置为 8 GB

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#默认情况下,cellchat使用object@data.signaling进行网络推断
#同时也提供了projectData函数,通过扩散过程基于高置信度实验验证的蛋白质互作网络中的邻近节点对基因表达值进行平滑处理。该功能在处理测序深度较浅的单细胞数据时尤为有用，因其能减少信号基因（特别是配体/受体亚基可能存在的零表达）的dropout效应。不担心其可能在扩散过程引入伪影，因其仅会引发极微弱的通讯信号。
# 原来是projectData，新版是smoothData函数
cellchat <- projectData(cellchat, adj = PPI.human)
# 假设你使用 meta$labels 创建了 cellchat


# 重建 cellchat 对象或更新 idents

options(future.globals.maxSize = 9 * 1024^3)  # 4 GB
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)
library(qs)
qsave(cellchat,"/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_middlecellchat.qs")
tumor_middle  <-qread("/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_middlecellchat.qs")
df.net.middle <- read_csv('/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_middledf.net.csv')
df.net.middle <- df.net.middle[,-1]
save(df.net,file = "/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_middledf.net.Rdata")
write.csv(df.net,"/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_middledf.net.csv")
# 计算聚合细胞-细胞通信网络
# 互作网络整合,可以设置soure和target，不设置就是默认全部
cellchat <- aggregateNet(cellchat)

# 可视化
groupSize <- as.numeric(table(cellchat@idents)) 
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFtumor_middleCellChat_netVisual1.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFtumor_middleCellChat_netVisual2.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
count_matrix <- cellchat@net[["count"]]

# 使用 pheatmap 绘制过滤后的矩阵
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFtumor_middleCellChat_heatmap.pdf", width = 8, height = 8)

pheatmap::pheatmap(count_matrix, 
                   border_color = "black", 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   fontsize = 10,
                   display_numbers = TRUE, 
                   number_color = "black", 
                   number_format = "%.0f")
dev.off()


caf_types = as.vector(unique(meta$labels))
caf_types = caf_types[-1]
caf_types
# 设置路径和参数
output_dir <- "/home/data/sdzl14/NSCLC/zong/fig/stromal"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

top_n <- 10

# 封装绘图函数
source_type = 'Tumor'
target_type = 'PI16+ proCAF'
plot_top_pathways <- function(df.net, source_type, target_type, output_dir, top_n = 10) {
  df_subset <- df.net %>%
    filter(source == source_type & target == target_type)
  
  if (nrow(df_subset) == 0) {
    message(paste0("No data for source='", source_type, "' and target='", target_type, "'. Skipping."))
    return()
  }
  
  pathway_weights <- df_subset %>%
    group_by(pathway_name) %>%
    summarise(total_prob = sum(prob), .groups = 'drop') %>%
    arrange(desc(total_prob))
  
  pathway_weights_topN <- pathway_weights %>%
    slice_max(order_by = total_prob, n = top_n) %>%
    arrange(desc(total_prob))
  
  p <- ggplot(pathway_weights_topN, aes(x = pathway_name, y = total_prob, fill = pathway_name)) +
    geom_bar(stat = "identity", color = "black") +
    coord_flip() +
    xlab("Pathway") +
    ylab("Total Interaction Weight (sum of prob)") +
    ggtitle(paste0("Top ", top_n, " Pathways in ", source_type, " -> ", target_type)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    ) +
    scale_fill_npg() +
    geom_text(aes(label = round(total_prob, 3)), vjust = -0.3, size = 3.5)
  
  filename <- file.path(output_dir, paste0("CellChat_tumor_middle", source_type, "_", target_type, ".pdf"))
  ggsave(filename, plot = p, width = 12, height = 12)
  
  return(p)
}

# 情况1: Tumor -> CAF
for (caf in caf_types) {
  plot_top_pathways(df.net, source_type = "Tumor", target_type = caf, output_dir = output_dir, top_n = top_n)
}

# 情况2: CAF -> Tumor
for (caf in caf_types) {
  plot_top_pathways(df.net, source_type = caf, target_type = "Tumor", output_dir = output_dir, top_n = top_n)
}
subset <- subset(merged_object,ident = 'tumor_edge')
data.input <- GetAssayData(subset, slot = 'data') # normalized data matrix
meta <- subset@meta.data[,c("Sample","celltype_coarse","Origin")]
colnames(meta) <-  c("Sample","labels","Origin")
identical(rownames(meta),colnames(data.input))
# 构建cellchat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
levels(cellchat@idents)
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
cellchat@DB <- CellChatDB.use
# Only uses the Secreted Signaling from CellChatDB v1
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole databa
future::plan("multisession", workers = 10) # do parallel
cellchat@data.signaling[cellchat@data.signaling <= 0] <- 0.01  # 加一个 pseudo-count 避免 log(0)
options(future.globals.maxSize = 8 * 1024^3)  # 设置为 8 GB

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#默认情况下,cellchat使用object@data.signaling进行网络推断
#同时也提供了projectData函数,通过扩散过程基于高置信度实验验证的蛋白质互作网络中的邻近节点对基因表达值进行平滑处理。该功能在处理测序深度较浅的单细胞数据时尤为有用，因其能减少信号基因（特别是配体/受体亚基可能存在的零表达）的dropout效应。不担心其可能在扩散过程引入伪影，因其仅会引发极微弱的通讯信号。
# 原来是projectData，新版是smoothData函数
cellchat <- projectData(cellchat, adj = PPI.human)
# 假设你使用 meta$labels 创建了 cellchat


# 重建 cellchat 对象或更新 idents

options(future.globals.maxSize = 9 * 1024^3)  # 4 GB
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)
library(qs)
qsave(cellchat,"/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_edgecellchat.qs")
save(df.net,file = "/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_edgedf.net.Rdata")
write.csv(df.net,"/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_edgedf.net.csv")
# 计算聚合细胞-细胞通信网络
# 互作网络整合,可以设置soure和target，不设置就是默认全部
cellchat <- normal_adjacent
cellchat <- aggregateNet(cellchat)

# 可视化
groupSize <- as.numeric(table(cellchat@idents)) 
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFnormal_adjacentCellChat_netVisual1.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFnormal_adjacentCellChat_netVisual2.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
count_matrix <- cellchat@net[["count"]]

# 使用 pheatmap 绘制过滤后的矩阵
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFnormal_adjacentCellChat_heatmap.pdf", width = 8, height = 8)

pheatmap::pheatmap(count_matrix, 
                   border_color = "black", 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   fontsize = 10,
                   display_numbers = TRUE, 
                   number_color = "black", 
                   number_format = "%.0f")
dev.off()


caf_types = as.vector(unique(cellchat@meta[["labels"]]))
caf_types = caf_types[-1]
caf_types
# 设置路径和参数
output_dir <- "/home/data/sdzl14/NSCLC/zong/fig/stromal"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

top_n <- 10
df.net <- df.net.adjacent
library(ggsci)
# 筛选 source=Tumor, target=Macro 的行
df_macro <- df.net %>%
  filter(source == "ADH1B+ Alveolar_CAF" & target == "Tumor")

# 按 pathway_name 分组，汇总 prob 权重
pathway_weights <- df_macro %>%
  group_by(pathway_name) %>%
  summarise(total_prob = sum(prob)) %>%
  arrange(desc(total_prob))
top_n <- 10
pathway_weights_topN <- pathway_weights %>%
  slice_max(order_by = total_prob, n = top_n) %>%
  arrange(desc(total_prob))

library(ggplot2)
library(ggthemes)
library(ggsci)

# 设置颜色方案
my_colors <- brewer.pal(n = 9, name = "Set1")[1:length(unique(pathway_weights_topN$pathway_name))]

# 横向柱状图，更好看也更清晰
ggplot(pathway_weights_topN, aes(x = pathway_name, y = total_prob, fill = pathway_name)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +  # 横向显示
  xlab("Pathway") +
  ylab("Total Interaction Weight (sum of prob)") +
  ggtitle(paste0("Top ", top_n, " Pathways in ADH1B+ Alveolar_CAF -> Tumor")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"  # 不显示图例，因为颜色只为美观
  ) +
  scale_fill_npg()  +
  geom_text(aes(label = round(total_prob, 3)), vjust = -0.3, size = 3.5)
ggsave(filename = "/home/data/sdzl14/NSCLC/zong/fig/stromal/CellChat_normal_adjacent_ADH1BAlveolar_CAF_Tumor.pdf",width = 12,height = 12)
levels(cellchat@idents) 
cellchat <- tumor_metastasis
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFMACRO/CellChat_tumor_metastasis_SEMA3_Tumor_netVisual_bubble.pdf",width = 9,height = 4.5)
netVisual_bubble(cellchat, sources  = c("ADH1B+ Alveolar_CAF" ,"COL11A1+ myCAF" ,"COL15A1+ mCAF", "LAMA2+ myCAF","MSLN+ iCAF" ,"PI16+ proCAF","POSTN+ mCAF"), 
                 targets.use  = c("Macro_CCL18" , "Macro_CHI3L1",'Macro_CXCL3',"Macro_CXCL9","Macro_SELENOP","Macro_SPP1" ), 
                 signaling = c("SEMA3"),
                 remove.isolate = T)
dev.off()
ggsave("/home/data/sdzl14/NSCLC/zong/fig/stromal/CellChat_tumor_edge_POSTNmCAF_Tumor_netVisual_bubble.pdf",width = 9,height = 6)

subset <- subset(merged_object,ident = 'normal_adjacent')
data.input <- GetAssayData(subset, slot = 'data') # normalized data matrix
meta <- subset@meta.data[,c("Sample","celltype_coarse","Origin")]
colnames(meta) <-  c("Sample","labels","Origin")
identical(rownames(meta),colnames(data.input))
# 构建cellchat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
levels(cellchat@idents)
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
cellchat@DB <- CellChatDB.use
# Only uses the Secreted Signaling from CellChatDB v1
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole databa
future::plan("multisession", workers = 10) # do parallel
cellchat@data.signaling[cellchat@data.signaling <= 0] <- 0.01  # 加一个 pseudo-count 避免 log(0)
options(future.globals.maxSize = 8 * 1024^3)  # 设置为 8 GB

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#默认情况下,cellchat使用object@data.signaling进行网络推断
#同时也提供了projectData函数,通过扩散过程基于高置信度实验验证的蛋白质互作网络中的邻近节点对基因表达值进行平滑处理。该功能在处理测序深度较浅的单细胞数据时尤为有用，因其能减少信号基因（特别是配体/受体亚基可能存在的零表达）的dropout效应。不担心其可能在扩散过程引入伪影，因其仅会引发极微弱的通讯信号。
# 原来是projectData，新版是smoothData函数
cellchat <- projectData(cellchat, adj = PPI.human)
# 假设你使用 meta$labels 创建了 cellchat


# 重建 cellchat 对象或更新 idents

options(future.globals.maxSize = 9 * 1024^3)  # 4 GB
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)
library(qs)
qsave(cellchat,"/home/data/sdzl14/NSCLC/zong/stromal/CAFnormal_adjacentcellchat.qs")
save(df.net,file = "/home/data/sdzl14/NSCLC/zong/stromal/CAFnormal_adjacentdf.net.Rdata")
write.csv(df.net,"/home/data/sdzl14/NSCLC/zong/stromal/CAFnormal_adjacentdf.net.csv")
# 计算聚合细胞-细胞通信网络
# 互作网络整合,可以设置soure和target，不设置就是默认全部
cellchat <- aggregateNet(cellchat)

# 可视化
groupSize <- as.numeric(table(cellchat@idents)) 
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFnormal_adjacentCellChat_netVisual1.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFnormal_adjacentCellChat_netVisual2.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
count_matrix <- cellchat@net[["count"]]

# 使用 pheatmap 绘制过滤后的矩阵
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFnormal_adjacentCellChat_heatmap.pdf", width = 8, height = 8)

pheatmap::pheatmap(count_matrix, 
                   border_color = "black", 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   fontsize = 10,
                   display_numbers = TRUE, 
                   number_color = "black", 
                   number_format = "%.0f")
dev.off()


caf_types = as.vector(unique(meta$labels))
caf_types = caf_types[-1]
caf_types
# 设置路径和参数
output_dir <- "/home/data/sdzl14/NSCLC/zong/fig/stromal"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

top_n <- 10

# 封装绘图函数
plot_top_pathways <- function(df.net, source_type, target_type, output_dir, top_n = 10) {
  df_subset <- df.net %>%
    filter(source == source_type & target == target_type)
  
  if (nrow(df_subset) == 0) {
    message(paste0("No data for source='", source_type, "' and target='", target_type, "'. Skipping."))
    return()
  }
  
  pathway_weights <- df_subset %>%
    group_by(pathway_name) %>%
    summarise(total_prob = sum(prob), .groups = 'drop') %>%
    arrange(desc(total_prob))
  
  pathway_weights_topN <- pathway_weights %>%
    slice_max(order_by = total_prob, n = top_n) %>%
    arrange(desc(total_prob))
  
  p <- ggplot(pathway_weights_topN, aes(x = pathway_name, y = total_prob, fill = pathway_name)) +
    geom_bar(stat = "identity", color = "black") +
    coord_flip() +
    xlab("Pathway") +
    ylab("Total Interaction Weight (sum of prob)") +
    ggtitle(paste0("Top ", top_n, " Pathways in ", source_type, " -> ", target_type)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    ) +
    scale_fill_npg() +
    geom_text(aes(label = round(total_prob, 3)), vjust = -0.3, size = 3.5)
  
  filename <- file.path(output_dir, paste0("CellChat", source_type, "_", target_type, ".pdf"))
  ggsave(filename, plot = p, width = 12, height = 12)
  
  return(p)
}

# 情况1: Tumor -> CAF
for (caf in caf_types) {
  plot_top_pathways(df.net, source_type = "Tumor", target_type = caf, output_dir = output_dir, top_n = top_n)
}

# 情况2: CAF -> Tumor
for (caf in caf_types) {
  plot_top_pathways(df.net, source_type = caf, target_type = "Tumor", output_dir = output_dir, top_n = top_n)
}
subset <- subset(merged_object,ident = 'tumor_metastasis')
data.input <- GetAssayData(subset, slot = 'data') # normalized data matrix
meta <- subset@meta.data[,c("Sample","celltype_coarse","Origin")]
colnames(meta) <-  c("Sample","labels","Origin")
identical(rownames(meta),colnames(data.input))
# 构建cellchat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
levels(cellchat@idents)
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
cellchat@DB <- CellChatDB.use
# Only uses the Secreted Signaling from CellChatDB v1
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole databa
future::plan("multisession", workers = 10) # do parallel
cellchat@data.signaling[cellchat@data.signaling <= 0] <- 0.01  # 加一个 pseudo-count 避免 log(0)
options(future.globals.maxSize = 8 * 1024^3)  # 设置为 8 GB

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#默认情况下,cellchat使用object@data.signaling进行网络推断
#同时也提供了projectData函数,通过扩散过程基于高置信度实验验证的蛋白质互作网络中的邻近节点对基因表达值进行平滑处理。该功能在处理测序深度较浅的单细胞数据时尤为有用，因其能减少信号基因（特别是配体/受体亚基可能存在的零表达）的dropout效应。不担心其可能在扩散过程引入伪影，因其仅会引发极微弱的通讯信号。
# 原来是projectData，新版是smoothData函数
cellchat <- projectData(cellchat, adj = PPI.human)
# 假设你使用 meta$labels 创建了 cellchat


# 重建 cellchat 对象或更新 idents

options(future.globals.maxSize = 9 * 1024^3)  # 4 GB
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)
library(qs)
qsave(cellchat,"/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_metastasiscellchat.qs")
save(df.net,file = "/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_metastasisdf.net.Rdata")
write.csv(df.net,"/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_metastasisdf.net.csv")
# 计算聚合细胞-细胞通信网络
# 互作网络整合,可以设置soure和target，不设置就是默认全部
cellchat <- aggregateNet(cellchat)

# 可视化
groupSize <- as.numeric(table(cellchat@idents)) 
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFtumor_metastasisCellChat_netVisual1.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFtumor_metastasisCellChat_netVisual2.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
count_matrix <- cellchat@net[["count"]]

# 使用 pheatmap 绘制过滤后的矩阵
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFtumor_metastasisCellChat_heatmap.pdf", width = 8, height = 8)

pheatmap::pheatmap(count_matrix, 
                   border_color = "black", 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   fontsize = 10,
                   display_numbers = TRUE, 
                   number_color = "black", 
                   number_format = "%.0f")
dev.off()


caf_types = as.vector(unique(meta$labels))
caf_types = caf_types[-1]
caf_types
# 设置路径和参数
output_dir <- "/home/data/sdzl14/NSCLC/zong/fig/stromal"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

top_n <- 10

# 封装绘图函数
plot_top_pathways <- function(df.net, source_type, target_type, output_dir, top_n = 10) {
  df_subset <- df.net %>%
    filter(source == source_type & target == target_type)
  
  if (nrow(df_subset) == 0) {
    message(paste0("No data for source='", source_type, "' and target='", target_type, "'. Skipping."))
    return()
  }
  
  pathway_weights <- df_subset %>%
    group_by(pathway_name) %>%
    summarise(total_prob = sum(prob), .groups = 'drop') %>%
    arrange(desc(total_prob))
  
  pathway_weights_topN <- pathway_weights %>%
    slice_max(order_by = total_prob, n = top_n) %>%
    arrange(desc(total_prob))
  
  p <- ggplot(pathway_weights_topN, aes(x = pathway_name, y = total_prob, fill = pathway_name)) +
    geom_bar(stat = "identity", color = "black") +
    coord_flip() +
    xlab("Pathway") +
    ylab("Total Interaction Weight (sum of prob)") +
    ggtitle(paste0("Top ", top_n, " Pathways in ", source_type, " -> ", target_type)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    ) +
    scale_fill_npg() +
    geom_text(aes(label = round(total_prob, 3)), vjust = -0.3, size = 3.5)
  
  filename <- file.path(output_dir, paste0("CellChat_tumor_metastasis", source_type, "_", target_type, ".pdf"))
  ggsave(filename, plot = p, width = 12, height = 12)
  
  return(p)
}

# 情况1: Tumor -> CAF
for (caf in caf_types) {
  plot_top_pathways(df.net, source_type = "Tumor", target_type = caf, output_dir = output_dir, top_n = top_n)
}

# 情况2: CAF -> Tumor
for (caf in caf_types) {
  plot_top_pathways(df.net, source_type = caf, target_type = "Tumor", output_dir = output_dir, top_n = top_n)
}
library(Seurat)
library(CellChat)
library(tidyverse)
library(patchwork)
library(scPRIT)
library(ggsci)
immune <- readRDS('/home/data/sdzl14/NSCLC/zong/immune.rds')
Idents(immune) <- 'immune_celltype_coarse'
macro <- subset(immune,idents = 'Macro')
macro@meta.data$Celltype <- macro@meta.data$immune_celltype
CAF <- readRDS('/home/data/sdzl14/NSCLC/zong/CAF.rds')
CAF@meta.data$Celltype <- CAF@meta.data$CAF_type 
merged_object <- merge(x = CAF, y = macro,  add.cell.ids = c("CAF", "macro"))
Idents(merged_object) <- 'Celltype'

Idents(merged_object) <- 'Tissue'
subset <- subset(merged_object,ident = 'tumor_metastasis')
NormalizeData(subset)
ScaleData(subset)
data.input <- GetAssayData(subset, slot = 'data') # normalized data matrix
meta <- subset@meta.data[,c("Sample","Celltype","Origin")]
colnames(meta) <-  c("Sample","labels","Origin")
identical(rownames(meta),colnames(data.input))
# 构建cellchat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
levels(cellchat@idents)
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
cellchat@DB <- CellChatDB.use
# Only uses the Secreted Signaling from CellChatDB v1
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole databa
future::plan("multisession", workers = 10) # do parallel
cellchat@data.signaling[cellchat@data.signaling <= 0] <- 0.01  # 加一个 pseudo-count 避免 log(0)
options(future.globals.maxSize = 8 * 1024^3)  # 设置为 8 GB

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#默认情况下,cellchat使用object@data.signaling进行网络推断
#同时也提供了projectData函数,通过扩散过程基于高置信度实验验证的蛋白质互作网络中的邻近节点对基因表达值进行平滑处理。该功能在处理测序深度较浅的单细胞数据时尤为有用，因其能减少信号基因（特别是配体/受体亚基可能存在的零表达）的dropout效应。不担心其可能在扩散过程引入伪影，因其仅会引发极微弱的通讯信号。
# 原来是projectData，新版是smoothData函数
cellchat <- projectData(cellchat, adj = PPI.human)
# 假设你使用 meta$labels 创建了 cellchat


# 重建 cellchat 对象或更新 idents

options(future.globals.maxSize = 9 * 1024^3)  # 4 GB
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)
library(qs)
qsave(cellchat,"/home/data/sdzl14/NSCLC/zong/stromal/CAFMACROtumor_metastasiscellchat.qs")
save(df.net,file = "/home/data/sdzl14/NSCLC/zong/stromal/CAFMACROtumor_metastasisdf.net.Rdata")
write.csv(df.net,"/home/data/sdzl14/NSCLC/zong/stromal/CAFMACROtumor_metastasisdf.net.csv")
tumor_middle  <-qread("/home/data/sdzl14/NSCLC/zong/stromal/CAFMACROtumor_middlecellchat.qs")
df.net.middle <- read_csv('/home/data/sdzl14/NSCLC/zong/stromal/CAFMACROtumor_middledf.net.csv')
df.net.middle <- df.net.middle[,-1]
tumor_edge  <-qread("/home/data/sdzl14/NSCLC/zong/stromal/CAFMACROtumor_edgecellchat.qs")
df.net.edge <- read_csv('/home/data/sdzl14/NSCLC/zong/stromal/CAFMACROtumor_edgedf.net.csv')
df.net.edge <- df.net.edge[,-1]
normal_adjacent  <-qread("/home/data/sdzl14/NSCLC/zong/stromal/CAFMACROnormal_adjacentcellchat.qs")
df.net.adjacent <- read_csv('/home/data/sdzl14/NSCLC/zong/stromal/CAFMACROnormal_adjacentdf.net.csv')
df.net.adjacent <- df.net.adjacent[,-1]
tumor_metastasis  <-qread("/home/data/sdzl14/NSCLC/zong/stromal/CAFMACROtumor_metastasiscellchat.qs")
df.net.metastasis <- read_csv('/home/data/sdzl14/NSCLC/zong/stromal/CAFMACROtumor_metastasisdf.net.csv')
df.net.metastasis <- df.net.metastasis[,-1]

# 提取每个 dataframe 中符合条件的 prob 总和
sum_metastasis <- sum(subset(df.net.metastasis,  pathway_name == "SPP1" &  target == 'Macro_CHI3L1')$prob)
sum_middle <- sum(subset(df.net.middle,  pathway_name == "SPP1"&  target == 'Macro_CHI3L1')$prob)
sum_edge <- sum(subset(df.net.edge, pathway_name == "SPP1" &  target == 'Macro_CHI3L1')$prob)
sum_adjacent <- sum(subset(df.net.adjacent,  pathway_name == "SPP1" &  target == 'Macro_CHI3L1')$prob)

# 构建绘图用的 dataframe
plot_data <- data.frame(
  group = c( "middle", "edge", "adjacent","metastasis"),
  total_prob = c( sum_middle, sum_edge, sum_adjacent,sum_metastasis)
)
# 手动指定每个柱子的颜色
plot_data$color <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")

# 设置 group 为因子，并指定顺序
plot_data$group <- factor(plot_data$group,
                          levels = c("middle", "edge", "adjacent", "metastasis"))
ggplot(plot_data, aes(x = group, y = total_prob, fill = group)) +
  
  # 使用自定义颜色映射
  geom_bar(stat = "identity", alpha = 0.8) +
  
  # 设置颜色
  scale_fill_manual(values = setNames(plot_data$color, plot_data$group)) +
  
  # 折线图连接各组
  geom_line(group = 1, color = "grey30", size = 1.2) +
  geom_point(size = 3, color = "grey30") +
  
  # 添加数值标签
  geom_text(aes(label = sprintf("%.4f", total_prob)), vjust = -0.5, size = 3.5) +
  
  # 标题与轴标签
  labs(
    title = "SPP1 Prob Tumor Sum Across Groups",
    x = "Group",
    y = "Total Prob"
  ) +
  
  # 主题美化
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"  # 隐藏图例
  )
ggsave("SPP1 targetMacroprob_comparison.pdf", width = 5, height = 5, dpi = 300)
# 定义函数：提取某个 df 中符合条件的 pathway prob 汇总
summarize_pathway <- function(df, group_name) {
  df %>%
    filter(target == "Macro_CHI3L1") %>%
    group_by(pathway_name) %>%
    summarise(total_prob = sum(prob, na.rm = TRUE)) %>%
    mutate(group = group_name)
}

# 分别处理四个 dataframe
pathway_metastasis <- summarize_pathway(df.net.metastasis, "metastasis")
pathway_middle     <- summarize_pathway(df.net.middle, "middle")
pathway_edge       <- summarize_pathway(df.net.edge, "edge")
pathway_adjacent   <- summarize_pathway(df.net.adjacent, "adjacent")

# 合并成一个完整的分析用 dataframe
pathway_comparison <- bind_rows(
  pathway_metastasis,
  pathway_middle,
  pathway_edge,
  pathway_adjacent)
library(tidyr)
library(ggplot2)

# 转换为宽格式
pathway_wide <- pathway_comparison %>%
  pivot_wider(names_from = group, values_from = total_prob)

# 绘制热图
head(df.net)[1:5,1:5]
ggplot(pathway_comparison, aes(x = group, y = pathway_name, fill = total_prob)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#F0E442", high = "#009E73", na.value = "lightgray") +
  theme_minimal() +
  labs(title = "TGFB Prob by Pathway Across Groups",
       x = "Group",
       y = "Pathway Name",
       fill = "Total Prob") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5))
# 定义函数：将每个 group 中 <10% 的合并为 Other
summarize_other <- function(df, group_col, value_col, label_col, threshold = 5) {
  df %>%
    group_by({{ group_col }}) %>%
    arrange(desc({{ value_col }}), .by_group = TRUE) %>%
    mutate(
      above_threshold = ifelse({{ value_col }} >= threshold, "Keep", "Other"),
      label = ifelse(above_threshold == "Keep", {{ label_col }}, "Other")
    ) %>%
    group_by({{ group_col }}, label) %>%
    summarise(value = sum({{ value_col }}, na.rm = TRUE)) %>%
    ungroup()
}
library(dplyr)
library(tidyr)
library(janitor)
??row_to_names
# 先转置数据，方便按列计算比例
pathway_wide_t <- pathway_wide %>%
  t() %>%
  as.data.frame() %>%
  row_to_names(row_number = 1) %>%
  mutate(rowname = rownames(.)) %>%
  pivot_longer(cols = -rowname, names_to = "pathway_name", values_to = "total_prob") %>%
  pivot_wider(names_from = rowname, values_from = total_prob)

# 按列计算百分比
pathway_percent <- pathway_wide %>%
  pivot_longer(cols = -pathway_name, names_to = "group", values_to = "prob") %>%
  group_by(group) %>%
  mutate(total = sum(prob, na.rm = TRUE),
         percent = (prob / total) * 100) %>%
  ungroup()
# 应用函数
pathway_percent_filtered <- summarize_other(pathway_percent,
                                            group_col = group,
                                            value_col = percent,
                                            label_col = pathway_name,
                                            threshold = 10)
library(ggplot2)
library(scales)
library(RColorBrewer)

# 获取唯一的 pathway_name 标签用于配色
labels <- unique(pathway_percent_filtered$label)
num_colors <- length(labels)

# 使用 RColorBrewer 的 Set3 配色方案（适合分类）
colors <- brewer.pal(n = num_colors, name = "Set3")

# 如果你想要自定义颜色，也可以手动指定：
# colors <- c("TGFb" = "#66C2A5", "Wnt" = "#FC8D62", "Other" = "#E78AC3")
pathway_percent_filtered$group <- factor(pathway_percent_filtered$group,
                                         levels = c("middle", "edge", "adjacent", "metastasis"))

# 绘图
p <- ggplot(pathway_percent_filtered, aes(x = group, y = value, fill = label)) +
  
  # 堆叠柱状图
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  
  # 添加百分比标签
  geom_text(
    aes(label = percent(value / 100)),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  
  # 设置配色
  scale_fill_manual(values = setNames(colors, labels)) +
  
  # 坐标轴和标题
  labs(
    title = "Pathway Distribution by Group (Percent > 10%)",
    x = "Group",
    y = "Percentage (%)",
    fill = "Pathway"
  ) +
  
  # Y轴格式化为百分比
  scale_y_continuous(labels = percent_format()) +
  
  # 主题美化
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "gray95", color = "black")
  )

# 显示图形
print(p)

# 保存图像（可选）
ggsave("TO_Macro_CHI3L1-pathway_distribution_percentage_beautiful.pdf", p, width = 7, height = 7, dpi = 300)
# 提取 source 为 Tumor 的数据，并按 group 和 target 统计 prob 总和
df.net.metastasis <- df.net.metastasis[df.net.metastasis$source != "AM", ]
df.net.metastasis <- df.net.metastasis[df.net.metastasis$target != "AM", ]
plot_data <- bind_rows(
  list(
    middle = df.net.middle,
    edge = df.net.edge,
    adjacent = df.net.adjacent,
    metastasis = df.net.metastasis
  ),
  .id = "group"
) %>%
  filter(source == "ADH1B+ Alveolar_CAF") %>%
  group_by(group, target) %>%
  summarise(total_prob = sum(prob, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("middle", "edge", "adjacent", "metastasis")))
plot_data_percent <- plot_data %>%
  group_by(group) %>%
  mutate(total = sum(total_prob),
         percent = (total_prob / total) * 100) %>%
  ungroup()
plot_data_filtered <- plot_data_percent %>%
  group_by(group) %>%
  arrange(desc(percent), .by_group = TRUE) %>%
  mutate(label = ifelse(percent < 10, "Other", target)) %>%
  group_by(group, label) %>%
  summarise(percent = sum(percent)) %>%
  ungroup()
# 获取唯一 target 标签用于配色
labels <- unique(plot_data_filtered$label)
num_colors <- length(labels)

# 使用 RColorBrewer 配色方案
colors <- brewer.pal(n = num_colors, name = "Set3")

# 绘图
p <- ggplot(plot_data_filtered, aes(x = group, y = percent, fill = label)) +
  
  # 堆叠柱状图
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  
  # 添加百分比标签
  geom_text(
    aes(label = percent(percent / 100)),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  
  # 设置配色
  scale_fill_manual(values = setNames(colors, labels)) +
  
  # 坐标轴与标题
  labs(
    title = "Target Distribution by Group (source = ADH1B+ Alveolar_CAF, Percent > 10%)",
    x = "Group",
    y = "Percentage (%)",
    fill = "Target"
  ) +
  
  # Y轴格式化为百分比
  scale_y_continuous(labels = percent_format()) +
  
  # 主题美化
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "gray95", color = "black")
  )

# 显示图形
print(p)

# 保存图像
ggsave("ADH1BAlveolar_CAFto_distribution_percentage.pdf", p, width = 10, height = 6, dpi = 300)
library(dplyr)
library(tidyr)

# 合并四个组的 dataframe 并筛选 target == "Tumor"
plot_data <- bind_rows(
  list(
    middle = df.net.middle,
    edge = df.net.edge,
    adjacent = df.net.adjacent,
    metastasis = df.net.metastasis
  ),
  .id = "group"
) %>%
  filter(source == "Macro_CHI3L1") %>%  # 修改为 target == "Tumor"
  group_by(group, source) %>%
  summarise(total_prob = sum(prob, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("middle", "edge", "adjacent", "metastasis")))
plot_data_percent <- plot_data %>%
  group_by(group) %>%
  mutate(total = sum(total_prob),
         percent = (total_prob / total) * 100) %>%
  ungroup()
plot_data_filtered <- plot_data_percent %>%
  group_by(group) %>%
  arrange(desc(percent), .by_group = TRUE) %>%
  mutate(label = ifelse(percent < 2, "Other", source)) %>%
  group_by(group, label) %>%
  summarise(percent = sum(percent)) %>%
  ungroup()
# 获取唯一 source 标签用于配色
labels <- unique(plot_data_filtered$label)
num_colors <- length(labels)

# 使用 RColorBrewer 配色方案
colors <- brewer.pal(n = num_colors, name = "Set3")

# 绘图
p <- ggplot(plot_data_filtered, aes(x = group, y = percent, fill = label)) +
  
  # 堆叠柱状图
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  
  # 添加百分比标签
  geom_text(
    aes(label = percent(percent / 100)),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  
  # 设置配色
  scale_fill_manual(values = setNames(colors, labels)) +
  
  # 坐标轴与标题
  labs(
    title = "Source Distribution by Group (Target = Tumor, Percent > 10%)",
    x = "Group",
    y = "Percentage (%)",
    fill = "Source"
  ) +
  
  # Y轴格式化为百分比
  scale_y_continuous(labels = percent_format()) +
  
  # 主题美化
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "gray95", color = "black")
  )

# 显示图形
print(p)
dev.off()
# 保存图像
pdf("source_to_tumor_distribution_target_tumor.pdf", width = 7, height = 7)
ggsave("source_to_tumor_distribution_target_tumor.png", p, width = 10, height = 6, dpi = 300)
table(merged_object@meta.data[["immune_celltype"]])
Idents(merged_object) <- 'immune_celltype_coarse'
subset <- subset(merged_object , ident = c('Tumor','Macro'))
Idents(subset) <- 'Tissue'
tumor_middle <- subset(subset, ident = 'tumor_middle')
tumor_edge <- subset(subset, ident = 'tumor_edge')
tumor_metastasis <- subset(subset,ident = 'tumor_metastasis')
options(future.globals.maxSize = 4 * 1024^3) 
NormalizeData(tumor_metastasis)
ScaleData(tumor_metastasis)

data.input <- GetAssayData(tumor_metastasis, slot = 'data') # normalized data matrix
meta <- tumor_metastasis@meta.data[,c("Sample","immune_celltype","Origin")]



colnames(meta) <-  c("Sample","labels","Origin")
table(meta$labels)

## 根据研究情况进行细胞排序
meta$labels <- factor(meta$labels ,levels = celltype_order)
table(meta$labels)

# 根据 meta$labels 的顺序进行排序
ordered_indices <- order(meta$labels)
# 对 meta 和 data.input 进行排序
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta),colnames(data.input))



dim(cellchat)
head(cellchat@data.signaling)[1:5,1:5]
cellchat@data.signaling <- data.raw
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE)

# 构建cellchat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
levels(cellchat@idents)
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
cellchat@DB <- CellChatDB.use
# Only uses the Secreted Signaling from CellChatDB v1
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole databa
future::plan("multisession", workers = 10) # do parallel
cellchat@data.signaling[cellchat@data.signaling <= 0] <- 0.01  # 加一个 pseudo-count 避免 log(0)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#默认情况下,cellchat使用object@data.signaling进行网络推断
#同时也提供了projectData函数,通过扩散过程基于高置信度实验验证的蛋白质互作网络中的邻近节点对基因表达值进行平滑处理。该功能在处理测序深度较浅的单细胞数据时尤为有用，因其能减少信号基因（特别是配体/受体亚基可能存在的零表达）的dropout效应。不担心其可能在扩散过程引入伪影，因其仅会引发极微弱的通讯信号。
# 原来是projectData，新版是smoothData函数
cellchat <- projectData(cellchat, adj = PPI.human)
# 假设你使用 meta$labels 创建了 cellchat
levels(meta$labels)
unique(meta$labels)

# 清理未使用的 factor levels
meta$labels <- droplevels(meta$labels)

# 重建 cellchat 对象或更新 idents
cellchat@idents <- meta$labels

cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 
cellchat <- filterCommunication(cellchat, min.cells = 8)
cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)
library(qs)
qsave(cellchat,"/home/data/sdzl14/NSCLC/zong/immune/tumor_metastasis-macro.qs")
save(df.net,file = "/home/data/sdzl14/NSCLC/zong/immune/tumor_metastasis-macrodf.net.Rdata")
write.csv(df.net,"/home/data/sdzl14/NSCLC/zong/immune/tumor_metastasis-macrodf.net.csv")
d

