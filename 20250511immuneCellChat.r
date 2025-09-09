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
devtools::install_local("/Users/zhangjinyang/Desktop/CellChat-main.zip")#后面这个就是包的路径
devtools::install_local("~/NSCLC/CellChat-main-2.zip")
library(svglite)
data.input <- GetAssayData(epi_counts, slot = 'data') # normalized data matrix
meta <- scRNA@meta.data[,c("orig.ident","celltype")]
colnames(meta) <-  c("group","labels")
table(meta$labels)
meta$labels <- gsub(" cells", "", meta$labels)
#meta$labels <- sub("\\(.*\\)", "", meta$labels)
table(meta$labels)

identical(rownames(meta),colnames(data.input))
## 根据研究情况进行细胞排序
celltype_order <- c(
  "T/NK", 
  "Th1", 
  "Th17", 
  "Tm", 
  "Treg", 
  "Naive T", 
  "ELK4+T", 
  "ZNF793+T", 
  "ZSCAN12+T",
  "B", 
  "VSMCs", 
  "endothelial", 
  "epithelial/cancer", 
  "fibroblasts", 
  "mast", 
  "myeloid", 
  "plasma", 
  "proliferative"
)
meta$labels <- factor(meta$labels ,levels = celltype_order)
table(meta$labels)

# 根据 meta$labels 的顺序进行排序
ordered_indices <- order(meta$labels)
# 对 meta 和 data.input 进行排序
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta),colnames(data.input))

# 构建cellchat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
levels(cellchat@idents)
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
cellchat@DB <- CellChatDB.use
# Only uses the Secreted Signaling from CellChatDB v1
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole databa
future::plan("multisession", workers = 1) # do parallel
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
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)

qsave(cellchat,"cellchat.qs")
save(df.net,file = "df.net.Rdata")
write.csv(df.net,"df.net.csv")
# 互作网络整合,可以设置soure和target，不设置就是默认全部
cellchat <- aggregateNet(cellchat)
# 可视化
groupSize <- as.numeric(table(cellchat@idents)) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
pheatmap::pheatmap(cellchat@net$count, border_color = "black", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T,number_color="black",number_format = "%.0f")




anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:30)
tumor <- readRDS('/home/data/sdzl14/NSCLC/zong/malignant_integrated.rds')
tumor@meta.data$immune_celltype_coarse <- 'Tumor'
tumor@meta.data$immune_celltype <- 'Tumor'
immune <- readRDS('/home/data/sdzl14/NSCLC/zong/immune.rds')
dim(immune)
# 合并Seurat对象
merged_object <- merge(x = tumor, y = immune,  add.cell.ids = c("Tumor", "Immune"))
dim(merged_object)
idents
celltype_order <- c(
  "CD4", 
  "CD8", 
  "Cycling T", 
  "Treg", 
  "NK", 
  "cDC1", 
  "cDC2", 
  "mregDC",
  "pDC",
  "Cycling B", 
  "GC B",
  "Memory B", 
  "Naive B",
  "Plasma", 
  "Macro", 
  "Mono", 
  "Mast",
  "AM", 
  "Neutrophil", 
  "ILC",
  "Tumor"
)

Idents(merged_object) <- "immune_celltype_coarse"
table(Idents(merged_object))
merged_object <- subset(merged_object, ident = celltype_order)
Idents(merged_object) <- "Tissue"
levels(merged_object)
tumor_middle <- subset(merged_object, idents = "tumor_middle")
tumor_edge <- subset(merged_object, idents = "tumor_edge")
normal_adjacent <- subset(merged_object, idents = "normal_adjacent")
tumor_metastasis  <- subset(merged_object, idents = "tumor_metastasis")
options(future.globals.maxSize = 4 * 1024^3) 
NormalizeData(tumor_middle)
ScaleData(tumor_middle)

data.input <- GetAssayData(tumor_middle, slot = 'data') # normalized data matrix
meta <- tumor_middle@meta.data[,c("Sample","immune_celltype","immune_celltype_coarse","Origin")]



colnames(meta) <-  c("Sample","immune_celltype","labels","Origin")
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
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)
library(qs)
qsave(cellchat,"/home/data/sdzl14/NSCLC/zong/immune/tumor_middlecellchat.qs")
save(df.net,file = "/home/data/sdzl14/NSCLC/zong/immune/tumor_middledf.net.Rdata")
write.csv(df.net,"/home/data/sdzl14/NSCLC/zong/immune/tumor_middledf.net.csv")
 计算聚合细胞-细胞通信网络
# 互作网络整合,可以设置soure和target，不设置就是默认全部
 tumor_edge <- aggregateNet(tumor_edge)
tumor_edge <- qread('/home/data/sdzl14/NSCLC/zong/immune/tumor_edgecellchat.qs')
# 可视化
groupSize <- as.numeric(table(cellchat@idents)) 
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_middle/CellChat_netVisual1.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_middle/CellChat_netVisual2.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/CellChat_netVisual_circle.pdf", width = 8, height = 6,dpi =300)
dev.off()
excluded_cells <- c("Cycling T", "cDC1", "mregDC", "pDC", "Cycling B", "Naive B", "Neutrophil")
count_matrix <- tumor_edge@net[["count"]]

# 删除被排除的细胞类型的行和列
count_filtered <- count_matrix[!(rownames(count_matrix) %in% excluded_cells), 
                               !(colnames(count_matrix) %in% excluded_cells)]

# 使用 pheatmap 绘制过滤后的矩阵
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_edge/CellChat_heatmap.pdf", width = 8, height = 8)

pheatmap::pheatmap(count_filtered, 
                   border_color = "black", 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   fontsize = 10,
                   display_numbers = TRUE, 
                   number_color = "black", 
                   number_format = "%.0f")
dev.off()
# 生成颜色向量（例如使用彩虹色）
color.use <- rainbow(nrow(df.net))
# 将颜色向量命名为矩阵的行名
names(color.use) <- rownames(df.net)
df.net <- as.data.frame(cellchat@net$weight)
celltype_order <- c(

  "CD4", 
  "CD8", 
  "Cycling T", 
  "Treg", 
  "NK", 
  "cDC1", 
  "cDC2", 
  "mregDC",
  "pDC",
  "GC B",
  "Memory B", 
  "Naive B",
  "Plasma", 
  "Macro", 
  "Mono", 
  "Mast",
  "AM", 
  "Neutrophil", 
  "Tumor"
)
df.net <- df.net[celltype_order,]#行排序
df.net <- df.net[,celltype_order] %>% as.matrix()
dim(df.net)
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_middle/CellChat_frow.pdf",width = 12,height = 12)
# 如果图片显示不全,需要考虑是不是重新设置mfrow
n_plots <- nrow(df.net)
n_col <- 4
n_row <- ceiling(n_plots / n_col)

par(mfrow = c(n_row, n_col), xpd = TRUE, mar = c(1,1,1,1))
# 生成颜色向量（例如使用彩虹色）
color.use <- rainbow(nrow(mat))
# 将颜色向量命名为矩阵的行名
names(color.use) <- rownames(mat)
mat <- as.data.frame(cellchat@net$weight)
celltype_order
mat <- mat[celltype_order,]#行排序
mat <- mat[,celltype_order] %>% as.matrix()
# 筛选 source=Tumor, target=Macro 的行
df_macro <- df.net %>%
  filter(source == "Tumor" & target == "Tumor")

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
  ggtitle(paste0("Top ", top_n, " Pathways in Tumore -> Tumor")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"  # 不显示图例，因为颜色只为美观
  ) +
  scale_fill_npg()  +
  geom_text(aes(label = round(total_prob, 3)), vjust = -0.3, size = 3.5)
ggsave(filename = "/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_middle/CellChat_Tumor_tumor.pdf",width = 12,height = 12)
# 批量绘图
head(df.net)
for (i in 1:ncol(df.net)) {
  
  # 提取单条通讯关系（source -> target）
  df.single <- df.net[i, , drop = FALSE]
  
  # 检查是否为方阵
  if (nrow(df.single) != ncol(df.single)) {
    warning(paste("Row and column numbers do not match for row:", i))
    next
  }
  
  # 绘图
  netVisual_circle(
    df.single,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    arrow.size = 0.05,
    arrow.width = 1,
    edge.weight.max = max(df.net),
    title.name = rownames(df.net)[i],
    color.use = color.use
  )
}

# 关闭设备
dev.off()
# 需要指定source和target
# sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(tumor_middle@idents) 
?netVisual_bubble
netVisual_bubble(cellchat, sources.use = seq(1:4), 
                 targets.use = c(18), remove.isolate = TRUE)
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_middle/CellChat_SPP1bubbleplot.pdf",width = 7,height = 4)

# 还可以增加signaling参数用于展示特定的配受体
cellchat@netP$pathways
netVisual_bubble(tumor_middle, targets.use = 'Tumor', 
                 sources.use  = c(1:20), 
                 signaling = c('SPP1'),
                 remove.isolate = FALSE)
ggsave("bubbleplot2.pdf",width = 5,height = 10)







cellchat@netP$pathways
levels(cellchat@idents) 
cellchat@netP$pathways
pathways.show <- "VEGF"
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") 
# 获取所有信号通路名称
all_pathways <- cellchat@netP$pathways

# 定义颜色（可选）
color.use <- rainbow(length(all_pathways))

# 遍历所有通路
for (pathway in all_pathways) {
  
  # 设置输出文件名，以通路命名
  filename <- paste0("/home/data/sdzl14/NSCLC/zong/fig/immune/CellChat_", pathway, ".pdf")
  
  # 打开 PDF 设备
  pdf(filename, width = 12, height = 12)
  
  # 绘图：展示该通路的信号发送、接收和中介作用
  netAnalysis_signalingRole_network(
    cellchat,
    signaling = pathway,
    width = 8,
    height = 8,
    font.size = 10
  )
  
  # 关闭当前 PDF 设备
  dev.off()
}
pathways.show <- "CXCL"
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/CellChat_Hierarchy plot.pdf",width = 12,height = 12)
# Hierarchy plot
# vertex.receiver定义层次图的左边细胞
vertex.receiver = seq(1:9) # a numeric vector
netVisual_aggregate(cellchat, signaling = pathways.show,
                    vertex.receiver = vertex.receiver,layout= "hierarchy")
                  # vertex.size = groupSize)  
# circle plot
netVisual_aggregate(cellchat, signaling = pathways.show,layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# 分组弦图
group.cellType <- c(rep("T/NK", 9), "B","VSMCs","endothelial","epithelial/cancer",
                    "fibroblasts","mast","myeloid","plasma","proliferative" )
levels(cellchat@idents)
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show,
                     group = group.cellType,
                     title.name = paste0(pathways.show, " signaling network"))

# heatmap
par(mfrow=c(1,1))
pathways.show = 'SPP1'
p <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_middle/CellChat_SPP1.pdf",width = 6,height = 6,plot = p)

# 需要指定source和target
# sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(cellchat@idents) 
netVisual_bubble(cellchat, sources.use = seq(1:20), 
                 targets.use = c(21), remove.isolate = FALSE)
ggsave("bubbleplot_nont.pdf",width = 7,height = 20)



# 自定义signaling输入展示-所有通路汇总之后
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("MIF","SPP1","MK"))
netVisual_bubble(cellchat, sources.use = c(1:18),
                 targets.use = c(18), 
                 pairLR.use = pairLR.use,
                 remove.isolate = TRUE)
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_middle/bubbleplot-LR.pdf",width = 5,height = 10)
NormalizeData(tumor_edge)
ScaleData(tumor_edge)

data.input <- GetAssayData(tumor_edge, slot = 'data') # normalized data matrix
meta <- tumor_edge@meta.data[,c("Sample","immune_celltype","immune_celltype_coarse","Origin")]



colnames(meta) <-  c("Sample","immune_celltype","labels","Origin")
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
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)
library(qs)
qsave(cellchat,"/home/data/sdzl14/NSCLC/zong/immune/tumor_edgecellchat.qs")
save(df.net,file = "/home/data/sdzl14/NSCLC/zong/immune/tumor_edgedf.net.Rdata")
write.csv(df.net,"/home/data/sdzl14/NSCLC/zong/immune/tumor_edgedf.net.csv")
计算聚合细胞-细胞通信网络
# 互作网络整合,可以设置soure和target，不设置就是默认全部
cellchat <- aggregateNet(cellchat)
# 可视化
groupSize <- as.numeric(table(cellchat@idents)) 
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_edge/macroCellChat_netVisual1.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_edge/macroCellChat_netVisual2.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/CellChat_netVisual_circle.pdf", width = 8, height = 6,dpi =300)
dev.off()
excluded_cells <- c("Cycling T", "cDC1", "mregDC", "pDC", "GC B", "Naive B",'Neutrophil',"Cycling B")
count_matrix <- cellchat@net$count

# 删除被排除的细胞类型的行和列
count_filtered <- count_matrix[!(rownames(count_matrix) %in% excluded_cells), 
                               !(colnames(count_matrix) %in% excluded_cells)]

# 使用 pheatmap 绘制过滤后的矩阵
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_edge/macroCellChat_heatmap.pdf", width = 8, height = 8)

pheatmap::pheatmap(count_matrix, 
                   border_color = "black", 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   fontsize = 10,
                   display_numbers = TRUE, 
                   number_color = "black", 
                   number_format = "%.0f")
# 生成颜色向量（例如使用彩虹色）
color.use <- rainbow(nrow(df.net))
# 将颜色向量命名为矩阵的行名
names(color.use) <- rownames(df.net)
df.net <- as.data.frame(cellchat@net$weight)
celltype_order <- c(
  
  "CD4", 
  "CD8", 
  "Cycling T", 
  "Treg", 
  "NK", 
  "cDC1", 
  "cDC2", 
  "mregDC",
  "pDC",
  "GC B",
  "Memory B", 
  "Naive B",
  "Plasma", 
  "Macro", 
  "Mono", 
  "Mast",
  "AM", 
  "Neutrophil", 
  "Tumor"
)
df.net <- df.net[celltype_order,]#行排序
df.net <- df.net[,celltype_order] %>% as.matrix()
dim(df.net)
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_edge/CellChat_frow.pdf",width = 12,height = 12)
# 如果图片显示不全,需要考虑是不是重新设置mfrow
n_plots <- nrow(df.net)
n_col <- 4
n_row <- ceiling(n_plots / n_col)

par(mfrow = c(n_row, n_col), xpd = TRUE, mar = c(1,1,1,1))
# 生成颜色向量（例如使用彩虹色）
color.use <- rainbow(nrow(mat))
# 将颜色向量命名为矩阵的行名
names(color.use) <- rownames(mat)
mat <- as.data.frame(cellchat@net$weight)
celltype_order
mat <- mat[celltype_order,]#行排序
mat <- mat[,celltype_order] %>% as.matrix()
# 筛选 source=Tumor, target=Macro 的行
df_macro <- df.net %>%
  filter(source == "Tumor" & target == "Tumor")

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
  ggtitle(paste0("Top ", top_n, " Pathways in Tumor -> Tumor")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"  # 不显示图例，因为颜色只为美观
  ) +
  scale_fill_npg()  +
  geom_text(aes(label = round(total_prob, 3)), vjust = -0.3, size = 3.5)
ggsave(filename = "/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_edge/CellChat_Tumor_Tumor.pdf",width = 12,height = 12)
# 批量绘图
head(df.net)
for (i in 1:ncol(df.net)) {
  
  # 提取单条通讯关系（source -> target）
  df.single <- df.net[i, , drop = FALSE]
  
  # 检查是否为方阵
  if (nrow(df.single) != ncol(df.single)) {
    warning(paste("Row and column numbers do not match for row:", i))
    next
  }
  
  # 绘图
  netVisual_circle(
    df.single,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    arrow.size = 0.05,
    arrow.width = 1,
    edge.weight.max = max(df.net),
    title.name = rownames(df.net)[i],
    color.use = color.use
  )
}

# 关闭设备
dev.off()
# 需要指定source和target
# sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(cellchat@idents) 
?netVisual_bubble
netVisual_bubble(cellchat, sources.use = seq(1:4), 
                 targets.use = c(18), remove.isolate = TRUE)
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_edge/CellChat_MIFANNEXINbubbleplot.pdf",width = 7,height = 6)

# 还可以增加signaling参数用于展示特定的配受体
cellchat@netP$pathways
netVisual_bubble(cellchat, targets.use  = seq(1:19), 
                 sources.use  = c(14), 
                 signaling = c("MIF",'ANNEXIN'),
                 remove.isolate = FALSE)
ggsave("bubbleplot2.pdf",width = 5,height = 10)







cellchat@netP$pathways
levels(cellchat@idents) 
cellchat@netP$pathways
pathways.show <- "VEGF"
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") 
# 获取所有信号通路名称
all_pathways <- cellchat@netP$pathways

# 定义颜色（可选）
color.use <- rainbow(length(all_pathways))

# 遍历所有通路
for (pathway in all_pathways) {
  
  # 设置输出文件名，以通路命名
  filename <- paste0("/home/data/sdzl14/NSCLC/zong/fig/immune/CellChat_", pathway, ".pdf")
  
  # 打开 PDF 设备
  pdf(filename, width = 12, height = 12)
  
  # 绘图：展示该通路的信号发送、接收和中介作用
  netAnalysis_signalingRole_network(
    cellchat,
    signaling = pathway,
    width = 8,
    height = 8,
    font.size = 10
  )
  
  # 关闭当前 PDF 设备
  dev.off()
}
pathways.show <- "CXCL"
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/CellChat_Hierarchy plot.pdf",width = 12,height = 12)
# Hierarchy plot
# vertex.receiver定义层次图的左边细胞
vertex.receiver = seq(1:9) # a numeric vector
netVisual_aggregate(cellchat, signaling = pathways.show,
                    vertex.receiver = vertex.receiver,layout= "hierarchy")
# vertex.size = groupSize)  
# circle plot
netVisual_aggregate(cellchat, signaling = pathways.show,layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# 分组弦图
group.cellType <- c(rep("T/NK", 9), "B","VSMCs","endothelial","epithelial/cancer",
                    "fibroblasts","mast","myeloid","plasma","proliferative" )
levels(cellchat@idents)
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show,
                     group = group.cellType,
                     title.name = paste0(pathways.show, " signaling network"))

# heatmap
par(mfrow=c(1,1))
pathways.show = 'SPP1'
p <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_edge/CellChat_SPP1.pdf",width = 6,height = 6,plot = p)

# 需要指定source和target
# sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(cellchat@idents) 
netVisual_bubble(cellchat, sources.use = seq(1:20), 
                 targets.use = c(21), remove.isolate = FALSE)
ggsave("bubbleplot_nont.pdf",width = 7,height = 20)



# 自定义signaling输入展示-所有通路汇总之后
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("MIF","SPP1","MK"))
netVisual_bubble(cellchat, sources.use = c(1:18),
                 targets.use = c(18), 
                 pairLR.use = pairLR.use,
                 remove.isolate = TRUE)
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_edge/bubbleplot-LR.pdf",width = 5,height = 10)
normal_adjacent <- subset(subset,ident = 'normal_adjacent')
tumor_metastasis <- subset(subset,ident = 'tumor_metastasis')
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
qsave(cellchat,"/home/data/sdzl14/NSCLC/zong/immune/tumor_metastasis-macrocellchat.qs")
save(df.net,file = "/home/data/sdzl14/NSCLC/zong/immune/tumor_metastasis-macrodf.net.Rdata")
write.csv(df.net,"/home/data/sdzl14/NSCLC/zong/immune/tumor_metastasis-macrodf.net.csv")

# 计算聚合细胞-细胞通信网络
# 互作网络整合,可以设置soure和target，不设置就是默认全部
cellchat <- aggregateNet(cellchat)
# 可视化
groupSize <- as.numeric(table(cellchat@idents)) 
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/normal_adjacent/macroCellChat_netVisual1.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/normal_adjacent/macroCellChat_netVisual2.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off()
excluded_cells <- c("Macro_SPP1", "Macro_CXCL3")
count_matrix <- cellchat@net$count

# 删除被排除的细胞类型的行和列
count_filtered <- count_matrix[!(rownames(count_matrix) %in% excluded_cells), 
                               !(colnames(count_matrix) %in% excluded_cells)]

# 使用 pheatmap 绘制过滤后的矩阵
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/normal_adjacent/macroCellChat_heatmap-.pdf", width = 8, height = 8)

pheatmap::pheatmap(count_matrix, 
                   border_color = "black", 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   fontsize = 10,
                   display_numbers = TRUE, 
                   number_color = "black", 
                   number_format = "%.0f")
dev.off()
# 生成颜色向量（例如使用彩虹色）
color.use <- rainbow(nrow(df.net))
# 将颜色向量命名为矩阵的行名
names(color.use) <- rownames(df.net)
df.net <- as.data.frame(cellchat@net$weight)
celltype_order <- c(
  
  "CD4", 
  "CD8", 
  "Cycling T", 
  "Treg", 
  "NK", 
  "cDC1", 
  "cDC2", 
  "mregDC",
  "pDC",
  "GC B",
  "Memory B", 
  "Naive B",
  "Plasma", 
  "Macro", 
  "Mono", 
  "Mast",
  "AM", 
  "Neutrophil", 
  "Tumor"
)
df.net <- df.net[celltype_order,]#行排序
df.net <- df.net[,celltype_order] %>% as.matrix()
dim(df.net)
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/normal_adjacent/CellChat_frow.pdf",width = 12,height = 12)
# 如果图片显示不全,需要考虑是不是重新设置mfrow
n_plots <- nrow(df.net)
n_col <- 4
n_row <- ceiling(n_plots / n_col)

par(mfrow = c(n_row, n_col), xpd = TRUE, mar = c(1,1,1,1))
# 生成颜色向量（例如使用彩虹色）
color.use <- rainbow(nrow(mat))
# 将颜色向量命名为矩阵的行名
names(color.use) <- rownames(mat)
mat <- as.data.frame(cellchat@net$weight)
celltype_order
mat <- mat[celltype_order,]#行排序
mat <- mat[,celltype_order] %>% as.matrix()
# 筛选 source=Tumor, target=Macro 的行
df_macro <- df.net %>%
  filter(source == "Tumor" & target == "Macro_CHI3L1")

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
  ggtitle(paste0("Top ", top_n, " Pathways in tumor -> Macro_CHI3L1")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"  # 不显示图例，因为颜色只为美观
  ) +
  scale_fill_npg()  +
  geom_text(aes(label = round(total_prob, 3)), vjust = -0.3, size = 3.5)
ggsave(filename = "/home/data/sdzl14/NSCLC/zong/fig/immune/normal_adjacent/Macro_CHI3L1CellChat_tumor-macro.pdf",width = 12,height = 12)
# 批量绘图
head(df.net)
for (i in 1:ncol(df.net)) {
  
  # 提取单条通讯关系（source -> target）
  df.single <- df.net[i, , drop = FALSE]
  
  # 检查是否为方阵
  if (nrow(df.single) != ncol(df.single)) {
    warning(paste("Row and column numbers do not match for row:", i))
    next
  }
  
  # 绘图
  netVisual_circle(
    df.single,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    arrow.size = 0.05,
    arrow.width = 1,
    edge.weight.max = max(df.net),
    title.name = rownames(df.net)[i],
    color.use = color.use
  )
}

# 关闭设备
dev.off()
# 需要指定source和target
# sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(cellchat@idents) 
?netVisual_bubble
netVisual_bubble(cellchat, sources.use = seq(1:4), 
                 targets.use = c(18), remove.isolate = TRUE)
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/normal_adjacent/CellChat_UGRP1bubbleplot.pdf",width = 7,height = 5)

# 还可以增加signaling参数用于展示特定的配受体
cellchat@netP$pathways
netVisual_bubble(cellchat, sources.use  = c(20), 
                 targets.use  = c(1:19), 
                 signaling = c("UGRP1",'GRN','GAS','ANNEXIN'),
                 remove.isolate = FALSE)
ggsave("bubbleplot2.pdf",width = 5,height = 10)







cellchat@netP$pathways
levels(cellchat@idents) 
cellchat@netP$pathways
pathways.show <- "VEGF"
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") 
# 获取所有信号通路名称
all_pathways <- cellchat@netP$pathways

# 定义颜色（可选）
color.use <- rainbow(length(all_pathways))

# 遍历所有通路
for (pathway in all_pathways) {
  
  # 设置输出文件名，以通路命名
  filename <- paste0("/home/data/sdzl14/NSCLC/zong/fig/immune/CellChat_", pathway, ".pdf")
  
  # 打开 PDF 设备
  pdf(filename, width = 12, height = 12)
  
  # 绘图：展示该通路的信号发送、接收和中介作用
  netAnalysis_signalingRole_network(
    cellchat,
    signaling = pathway,
    width = 8,
    height = 8,
    font.size = 10
  )
  
  # 关闭当前 PDF 设备
  dev.off()
}
pathways.show <- "CXCL"
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/CellChat_Hierarchy plot.pdf",width = 12,height = 12)
# Hierarchy plot
# vertex.receiver定义层次图的左边细胞
vertex.receiver = seq(1:9) # a numeric vector
netVisual_aggregate(cellchat, signaling = pathways.show,
                    vertex.receiver = vertex.receiver,layout= "hierarchy")
# vertex.size = groupSize)  
# circle plot
netVisual_aggregate(cellchat, signaling = pathways.show,layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# 分组弦图
group.cellType <- c(rep("T/NK", 9), "B","VSMCs","endothelial","epithelial/cancer",
                    "fibroblasts","mast","myeloid","plasma","proliferative" )
levels(cellchat@idents)
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show,
                     group = group.cellType,
                     title.name = paste0(pathways.show, " signaling network"))

# heatmap
par(mfrow=c(1,1))
pathways.show = 'SPP1'
p <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/normal_adjacent/CellChat_SPP1.pdf",width = 6,height = 6,plot = p)

# 需要指定source和target
# sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(cellchat@idents) 
netVisual_bubble(cellchat, sources.use = seq(1:20), 
                 targets.use = c(21), remove.isolate = FALSE)
ggsave("bubbleplot_nont.pdf",width = 7,height = 20)



# 自定义signaling输入展示-所有通路汇总之后
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("MIF","SPP1","MK"))
netVisual_bubble(cellchat, sources.use = c(1:18),
                 targets.use = c(18), 
                 pairLR.use = pairLR.use,
                 remove.isolate = TRUE)
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/normal_adjacent/bubbleplot-LR.pdf",width = 5,height = 10)
Idents(tumor_metastasis) <- 'immune_celltype_coarse'
unique(tumor_metastasis$immune_celltype_coarse)
DimPlot(tumor_metastasis)
tumor_metastasis <- subset(tumor_metastasis,idents != "NA")
NormalizeData(tumor_metastasis)
ScaleData(tumor_metastasis)

data.input <- GetAssayData(tumor_metastasis, slot = 'data') # normalized data matrix
meta <- tumor_metastasis@meta.data[,c("Sample","immune_celltype","immune_celltype_coarse","Origin")]



colnames(meta) <-  c("Sample","immune_celltype","labels","Origin")
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
#获取非 NA 的细胞索引
valid_cells <- !is.na(cellchat@idents)

# 更新 cellchat 对象中的 idents 和 data.signaling 等关键数据
cellchat@idents <- cellchat@idents[valid_cells]
cellchat@data.signaling <- cellchat@data.signaling[, valid_cells]
# 假设你的原始数据是 data.input，meta 是 meta
# 先过滤掉 meta 中为 NA 的行
meta_clean <- meta[!is.na(meta$labels), , drop = FALSE]

# 同步过滤 data.input 的列
data.input_clean <- data.input[, rownames(meta_clean)]
# 重新创建 cellchat 对象
cellchat <- createCellChat(object = data.input_clean, meta = meta_clean, group.by = "labels")
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 
cellchat <- filterCommunication(cellchat, min.cells = 8)
cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)
library(qs)
qsave(cellchat,"/home/data/sdzl14/NSCLC/zong/immune/tumor_metastasiscellchat.qs")
save(df.net,file = "/home/data/sdzl14/NSCLC/zong/immune/tumor_metastasisdf.net.Rdata")
write.csv(df.net,"/home/data/sdzl14/NSCLC/zong/immune/tumor_metastasisdf.net.csv")
# 计算聚合细胞-细胞通信网络
# 互作网络整合,可以设置soure和target，不设置就是默认全部
cellchat <- aggregateNet(cellchat)
# 可视化
groupSize <- as.numeric(table(cellchat@idents)) 
groupSize <- groupSize[-17]
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_edge/macroCellChat_netVisual3.pdf")
netVisual_circle(count_filtered, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_edge/macroCellChat_netVisual4.pdf")
weight <- cellchat@net$weight
weight <- weight[-17,-17]
netVisual_circle(weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/CellChat_netVisual_circle.pdf", width = 8, height = 6,dpi =300)
dev.off()
excluded_cells <- c('AM')
count_matrix <- cellchat@net$count

# 删除被排除的细胞类型的行和列
count_filtered <- count_matrix[!(rownames(count_matrix) %in% excluded_cells), 
                               !(colnames(count_matrix) %in% excluded_cells)]

# 使用 pheatmap 绘制过滤后的矩阵
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_metastasis/CellChat_heatmap.pdf", width = 8, height = 8)

pheatmap::pheatmap(count_filtered, 
                   border_color = "black", 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   fontsize = 10,
                   display_numbers = TRUE, 
                   number_color = "black", 
                   number_format = "%.0f")
# 生成颜色向量（例如使用彩虹色）
color.use <- rainbow(nrow(df.net))
# 将颜色向量命名为矩阵的行名
names(color.use) <- rownames(df.net)
df.net <- as.data.frame(cellchat@net$weight)
celltype_order <- c(
  
  "CD4", 
  "CD8", 
  "Cycling T", 
  "Treg", 
  "NK", 
  "cDC1", 
  "cDC2", 
  "mregDC",
  "pDC",
  "GC B",
  "Memory B", 
  "Naive B",
  "Plasma", 
  "Macro", 
  "Mono", 
  "Mast",
  "AM", 
  "Neutrophil", 
  "Tumor"
)
df.net <- df.net[celltype_order,]#行排序
df.net <- df.net[,celltype_order] %>% as.matrix()
dim(df.net)
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_metastasis/CellChat_frow.pdf",width = 12,height = 12)
# 如果图片显示不全,需要考虑是不是重新设置mfrow
n_plots <- nrow(df.net)
n_col <- 4
n_row <- ceiling(n_plots / n_col)

par(mfrow = c(n_row, n_col), xpd = TRUE, mar = c(1,1,1,1))
# 生成颜色向量（例如使用彩虹色）
color.use <- rainbow(nrow(mat))
# 将颜色向量命名为矩阵的行名
names(color.use) <- rownames(mat)
mat <- as.data.frame(cellchat@net$weight)
celltype_order
mat <- mat[celltype_order,]#行排序
mat <- mat[,celltype_order] %>% as.matrix()
# 筛选 source=Tumor, target=Macro 的行
df_macro <- df.net %>%
  filter(source == "Macro" & target == "Tumor")

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
  ggtitle(paste0("Top ", top_n, " Pathways in Macro -> Tumor")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"  # 不显示图例，因为颜色只为美观
  ) +
  scale_fill_npg()  +
  geom_text(aes(label = round(total_prob, 3)), vjust = -0.3, size = 3.5)
ggsave(filename = "/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_metastasis/CellChat_Macro-Tumor.pdf",width = 12,height = 12)
# 批量绘图
head(df.net)
for (i in 1:ncol(df.net)) {
  
  # 提取单条通讯关系（source -> target）
  df.single <- df.net[i, , drop = FALSE]
  
  # 检查是否为方阵
  if (nrow(df.single) != ncol(df.single)) {
    warning(paste("Row and column numbers do not match for row:", i))
    next
  }
  
  # 绘图
  netVisual_circle(
    df.single,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    arrow.size = 0.05,
    arrow.width = 1,
    edge.weight.max = max(df.net),
    title.name = rownames(df.net)[i],
    color.use = color.use
  )
}

# 关闭设备
dev.off()
# 需要指定source和target
# sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(cellchat@idents) 
?netVisual_bubble
netVisual_bubble(cellchat, sources.use = seq(1:4), 
                 targets.use = c(18), remove.isolate = TRUE)
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_metastasis/CellChat_SPP1bubbleplot.pdf",width = 7,height = 7)

# 还可以增加signaling参数用于展示特定的配受体
cellchat@netP$pathways
netVisual_bubble(cellchat, targets.use  = c(19), 
                 sources.use  = c(1:16,18,19), 
                 signaling = c('SPP1'),
                 remove.isolate = FALSE)
ggsave("bubbleplot2.pdf",width = 5,height = 10)







cellchat@netP$pathways
levels(cellchat@idents) 
cellchat@netP$pathways
pathways.show <- "VEGF"
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") 
# 获取所有信号通路名称
all_pathways <- cellchat@netP$pathways

# 定义颜色（可选）
color.use <- rainbow(length(all_pathways))

# 遍历所有通路
for (pathway in all_pathways) {
  
  # 设置输出文件名，以通路命名
  filename <- paste0("/home/data/sdzl14/NSCLC/zong/fig/immune/CellChat_", pathway, ".pdf")
  
  # 打开 PDF 设备
  pdf(filename, width = 12, height = 12)
  
  # 绘图：展示该通路的信号发送、接收和中介作用
  netAnalysis_signalingRole_network(
    cellchat,
    signaling = pathway,
    width = 8,
    height = 8,
    font.size = 10
  )
  
  # 关闭当前 PDF 设备
  dev.off()
}
pathways.show <- "CXCL"
pdf("/home/data/sdzl14/NSCLC/zong/fig/immune/CellChat_Hierarchy plot.pdf",width = 12,height = 12)
# Hierarchy plot
# vertex.receiver定义层次图的左边细胞
vertex.receiver = seq(1:9) # a numeric vector
netVisual_aggregate(cellchat, signaling = pathways.show,
                    vertex.receiver = vertex.receiver,layout= "hierarchy")
# vertex.size = groupSize)  
# circle plot
netVisual_aggregate(cellchat, signaling = pathways.show,layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# 分组弦图
group.cellType <- c(rep("T/NK", 9), "B","VSMCs","endothelial","epithelial/cancer",
                    "fibroblasts","mast","myeloid","plasma","proliferative" )
levels(cellchat@idents)
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show,
                     group = group.cellType,
                     title.name = paste0(pathways.show, " signaling network"))

# heatmap
par(mfrow=c(1,1))
pathways.show = 'SPP1'
p <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_metastasis/CellChat_SPP1.pdf",width = 6,height = 6,plot = p)

# 需要指定source和target
# sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(cellchat@idents) 
netVisual_bubble(cellchat, sources.use = seq(1:20), 
                 targets.use = c(21), remove.isolate = FALSE)
ggsave("bubbleplot_nont.pdf",width = 7,height = 20)



# 自定义signaling输入展示-所有通路汇总之后
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("MIF","SPP1","MK"))
netVisual_bubble(cellchat, sources.use = c(1:18),
                 targets.use = c(18), 
                 pairLR.use = pairLR.use,
                 remove.isolate = TRUE)
ggsave("/home/data/sdzl14/NSCLC/zong/fig/immune/tumor_metastasis/bubbleplot-LR.pdf",width = 5,height = 10)
df.net.metastasis <- df.net
head(df.net.metastasis)
df.net.middle <- read_csv('~/NSCLC/zong/immune/tumor_middledf.net.csv')
df.net.middle <- df.net.middle[,-1]
df.net.edge <- read_csv('~/NSCLC/zong/immune/tumor_edgedf.net.csv')
df.net.edge <- df.net.edge[,-1]
df.net.adjacent <- read_csv('~/NSCLC/zong/immune/normal_adjacentdf.net.csv')
df.net.adjacent <- df.net.adjacent[,-1]
# 提取每个 dataframe 中符合条件的 prob 总和
sum_metastasis <- sum(subset(df.net.metastasis,  pathway_name == "GALECTIN" &  target == 'Tumor')$prob)
sum_middle <- sum(subset(df.net.middle,  pathway_name == "GALECTIN"&  target == 'Tumor')$prob)
sum_edge <- sum(subset(df.net.edge, pathway_name == "GALECTIN" &  target == 'Tumor')$prob)
sum_adjacent <- sum(subset(df.net.adjacent,  pathway_name == "GALECTIN" &  target == 'Tumor')$prob)

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
    title = "GALECTIN Prob Tumor Sum Across Groups",
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
ggsave("GALECTIN_Tumor_prob_comparison.pdf", width = 5, height = 5, dpi = 300)
# 定义函数：提取某个 df 中符合条件的 pathway prob 汇总
summarize_pathway <- function(df, group_name) {
  df %>%
    filter(target == "Macro", source == "Tumor") %>%
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
  labs(title = "SPP1 Prob by Pathway Across Groups",
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
                                            threshold = 5)
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
ggsave("Tumor-Macro-pathway_distribution_percentage_beautiful.pdf", p, width = 7, height = 7, dpi = 300)
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
  filter(Target == "Tumor") %>%
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
  mutate(label = ifelse(percent < 8, "Other", target)) %>%
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
    title = "Target Distribution by Group (Target = Tumor, Percent > 10%)",
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
ggsave("tumor_to_target_distribution_percentage.png", p, width = 10, height = 6, dpi = 300)
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
  filter(target == "Tumor") %>%  # 修改为 target == "Tumor"
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
  mutate(label = ifelse(percent < 8, "Other", source)) %>%
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
