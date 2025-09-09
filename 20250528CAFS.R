library(Seurat)
.libPaths('~/R/x86_64-pc-linux-gnu-library/4.2')
library(CellChat)
library(tidyverse)
library(patchwork)
library(qs)
cellchat <- qread('/home/data/sdzl14/NSCLC/zong/stromal/CAFtumor_middlecellchat.qs')
cellchat <- aggregateNet(cellchat)
unique(cellchat@idents)
# 可视化
groupSize <- as.numeric(table(cellchat@idents)) 
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFtumor_middlecCellChat_netVisual1.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("/home/d ata/sdzl14/NSCLC/zong/fig/stromal/CAFtumor_middlecCellChat_netVisual2.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off()
excluded_cells <- c("PI16+ proCAF", "LAMA2+ myCAF", "COL15A1+ mCAF")
count_matrix <- cellchat@net[["count"]]

# 删除被排除的细胞类型的行和列
count_filtered <- count_matrix[!(rownames(count_matrix) %in% excluded_cells), 
                               !(colnames(count_matrix) %in% excluded_cells)]

# 使用 pheatmap 绘制过滤后的矩阵
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFtumor_middleCellChat_heatmap.pdf", width = 8, height = 8)

pheatmap::pheatmap(count_filtered, 
                   border_color = "black", 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   fontsize = 10,
                   display_numbers = TRUE, 
                   number_color = "black", 
                   number_format = "%.0f")
dev.off()

dim(df.net)
# 筛选 source=Tumor, target=Macro 的行
df.net <- subsetCommunication(cellchat)
df_macro <- df.net %>%
  filter(target == "Tumor" & source == "POSTN+ mCAF")

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
  ggtitle(paste0("Top ", top_n, " Pathways in POSTN+ mCAF ->Tumor ")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"  # 不显示图例，因为颜色只为美观
  ) +
  scale_fill_npg()  +
  geom_text(aes(label = round(total_prob, 3)), vjust = -0.3, size = 3.5)
ggsave(filename = "/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFtumor_middlecCellCha_POSTNmCAFt_Tumor_.pdf",width = 12,height = 12)
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
cellchat <- tumor_middle
cellchat <- aggregateNet(cellchat)

# 可视化
groupSize <- as.numeric(table(cellchat@idents)) 
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFMACRO/CAFtumor_metastasisCellChat_netVisual1.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFMACRO/CAFtumor_metastasisCellChat_netVisual2.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
count_matrix <- cellchat@net[["count"]]

# 使用 pheatmap 绘制过滤后的矩阵
pdf("/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFMACRO/CAFtumor_metastasisCellChat_heatmap.pdf", width = 8, height = 8)

pheatmap::pheatmap(count_matrix, 
                   border_color = "black", 
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE, 
                   fontsize = 10,
                   display_numbers = TRUE, 
                   number_color = "black", 
                   number_format = "%.0f")
dev.off()
# 筛选 source=Tumor, target=Macro 的行
df.net <- df.net.middle
df_macro <- df.net %>%
  filter(source == "POSTN+ mCAF" & target == "Macro_CHI3L1")

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
  ggtitle(paste0("Top ", top_n, " Pathways in POSTN+ mCAF -> Macro_CHI3L1")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"  # 不显示图例，因为颜色只为美观
  ) +
  scale_fill_npg()  +
  geom_text(aes(label = round(total_prob, 3)), vjust = -0.3, size = 3.5)
ggsave(filename = "/home/data/sdzl14/NSCLC/zong/fig/stromal/CAFMACRO/CellChat_Macro_CHI3L1_POSTNmCAF_tumor_middles.pdf",width = 12,height = 12)
# sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(cellchat@idents) 
cellchat <- tumor_edge
netVisual_bubble(cellchat, sources.use  = c(4,5,6,7,8), 
                 targets.use  = c(1,2,3,9,10), 
                 signaling = c('SPP1'),
                 remove.isolate = FALSE)
ggsave("tumor_edgeSPP1bubbleplot2.pdf",width = 7,height = 5)
