.libPaths('~/R/x86_64-pc-linux-gnu-library/4.2')
.libPaths()
options(bitmapType='cairo')
BiocManager::install('GeneOverlap')
devtools::install_github('smorabit/hdWGCNA', ref='dev')
devtools::install_local(path = '~/NSCLC/hdWGCNA-dev.zip')
library(hdWGCNA)
library(Seurat)
library(qs)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(magrittr)
library(hdWGCNA)
theme_set(theme_cowplot())
pdf("~/dimplot.pdf")
print(p)
dev.off()
malignant_integrated <- readRDS(file = '~/NSCLC/zong/LUSC.rds')
Idents(malignant_integrated) <- 'leiden'

DimPlot(malignant_integrated)
NormalizeData(malignant_integrated)
ScaleData(malignant_integrated)
str(malignant_integrated@assays[["RNA"]]@data)
dim(malignant_integrated@assays[["RNA"]]@data)
anyNA(malignant_integrated@assays[["RNA"]]@data)
FindVariableFeatures(malignant_integrated)
Idents(malignant_integrated) <- 'Pathtype'
# 随机抽取30000个细胞的索引
sampled_cells <- sample(seq_len(ncol(malignant_integrated)), size = 30000, replace = FALSE)

# 创建新的Seurat对象，仅包含采样的细胞
malignant_integrated <- subset(malignant_integrated, cells = sampled_cells)
?SetupForWGCNA
seurat_obj <- SetupForWGCNA(
  malignant_integrated,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Malignant" # the name of the hdWGCNA experiment
)
length(seurat_obj@misc$Malignant$wgcna_genes)
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Timepoint","Tissue","Origin","Pathtype","leiden"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Tissue' # set the Idents of the metacell seurat object
)
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = 'Pre',      
  group.by='Timepoint', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'counts' # using normalized data
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
#可视化
plot_list <- PlotSoftPowers(seurat_obj)

wrap_plots(plot_list)#这里显然计算所得我们的最佳软阈值是6
?ConstructNetwork
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power=5,
  setDatExpr=FALSE,
  corType = "pearson",
  networkType = "signed",
  TOMType = "signed",
  detectCutHeight = 0.995,
  minModuleSize = 25,
  mergeCutHeight = 0.1,
  overwrite_tom = TRUE,
  tom_outdir = "TOM", # 输出文件夹
  tom_name = 'MAC' # name of the topoligical overlap matrix written to disk
)
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power = 5,
  setDatExpr = FALSE,
  corType = "pearson",
  networkType = "signed",
  TOMType = "signed",
  deepSplit = 2,                 # 减少初始分割深度
  detectCutHeight = 0.998,       # 增加剪枝阈值，减少模块数量
  minModuleSize = 10,            # 允许更小的模块被保留
  mergeCutHeight = 0.3,          # 更积极地合并模块
  overwrite_tom = TRUE,
  tom_outdir = "TOM",
  tom_name = 'LUAD'
)
?PlotDendrogram
PlotDendrogram(seurat_obj, main='MAC hdWGCNA Dendrogram')
ggsave(filename = '~/NSCLC/zong/fig/MALIGNANT/WCGNA_LUADPlotDendrogram.pdf',height = 10,width = 10,dpi = 300)
seurat_obj <- ScaleData(seurat_obj,features=VariableFeatures(seurat_obj))
?ModuleEigengenes
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  scale.model.use="linear",
  assay = NULL,
  pc_dim = 1)
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'Timepoint', group_name = 'Pre',corFnc="bicor",corOptions="use='p'",harmonized=TRUE,assay=NULL,slot="data"
)
# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=T)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj,ncol = 4)

p
ggsave(filename =  '~/NSCLC/zong/fig/MALIGNANT/WCGNA_LUADPlotKMEs.pdf',height = 8,width = 12,dpi = 300)
TOM <- hdWGCNA::GetTOM(seurat_obj)
hdWGCNA::ModuleNetworkPlot(seurat_obj)
# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
#> head(modules[,1:6])
#gene_name module     color  kME_grey kME_LUAD1 kME_LUAD2
#HES4          HES4  LUAD1 turquoise 0.3565626 0.3601537 0.3348734
#ISG15        ISG15  LUAD2     brown 0.7167839 0.6851797 0.7433319
#AGRN          AGRN  LUAD3    salmon 0.3451876 0.3524225 0.3391006
#C1orf159  C1orf159  LUAD1 turquoise 0.3229881 0.3168387 0.3209049
#ACAP3        ACAP3  LUAD1 turquoise 0.2219663 0.2299583 0.2025830
#AURKAIP1  AURKAIP1  LUAD4      blue 0.8768879 0.8246393 0.8730361
#> head(hub_df)
#gene_name module       kME
#1    HNRNPC  LUAD1 0.8702770
#2     RPL37  LUAD1 0.8513113
#3     SRSF9  LUAD1 0.8479055
#4     PDIA3  LUAD1 0.8449099
#5     RAP2B  LUAD1 0.8425380
#6      RAC1  LUAD1 0.8399496
saveRDS(seurat_obj, file='~/NSCLC/zong/LUADhdWGCNA_object.rds')
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=4)
ggsave(filename =  '~/NSCLC/zong/fig/MALIGNANT/WCGNA_LUADModuleFeaturePlot.pdf',height = 12,width = 12,dpi = 300)
# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=4)
ggsave(filename =  '~/NSCLC/zong/fig/MALIGNANT/WCGNA_LUADModuleFeaturePlotscores.pdf',height = 12,width = 12,dpi = 300)
seurat_obj$leiden <- do.call(rbind, strsplit(as.character(seurat_obj$leiden), ' '))[,1]

ModuleRadarPlot(
  seurat_obj,
  group.by = 'leiden',
  axis.label.size=4,
  grid.label.size=4
)
ggsave(filename =  '~/NSCLC/zong/fig/MALIGNANT/WCGNA_LUADModuleRadarPlot.pdf',height = 12,width = 12,dpi = 300)
# plot module correlagram
ModuleCorrelogram(seurat_obj)
ggsave(filename =  '~/NSCLC/zong/fig/MALIGNANT/WCGNA_LUADModuleCorrelogram.pdf',height = 8,width = 8,dpi = 300)
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
# plot with Seurat's DotPlot function

p <- DotPlot(seurat_obj, features=mods, group.by = 'leiden')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p
ggsave()
# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
head(hub_df)
net <- seurat_obj@misc[["Malignant"]][["wgcna_net"]]
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
geneTree = net$dendrograms[[1]]
# 将标签转换为颜色
mergedColors = labels2colors(net$colors)
# 绘制树状图和模块颜色图
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Modulecolors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,addTextGuide = TRUE)
TOM <- GetTOM(seurat_obj)       
hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
head(MEs)[1:5,1:5]

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "LUAD"
)
modules<-GetModules(seurat_obj)
head(modules[,1:6])
print(levels(modules$module))
p <- PlotKMEs(seurat_obj, 
              ncol=5,
              n_hubs = 10, 
              text_size = 2,
              plot_widths = c(3, 2)) 
ggsave(filename = '~/NSCLC/zong/fig/MALIGNANT/WGCNA_plotKMEs-.pdf',plot = p,height = 2.5,width = 15,dpi = 300)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)
# 将标签转换为颜色
mergedColors = labels2colors(net$colors)
# 绘制树状图和模块颜色图
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Modulecolors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,addTextGuide = TRUE)
ggsave(filename = '~/NSCLC/zong/fig/MALIGNANT/WGCNA_树状图和模块颜色图.pdf',height = 8,width = 8,dpi = 300)
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
ggsave(filename = '~/NSCLC/zong/fig/MALIGNANT/WGCNA_模块特征基因聚类树热图图.pdf',height = 8,width = 8,dpi = 300)
pdf(file="12_dendrogram.pdf",width=8, height=8)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()
cytoDir="CytoscapeInput"#【构建文件夹】
dir.create(cytoDir)
# 检查 wgcna_net 的内容
wgcna_net <- seurat_obj@misc[["Malignant"]][["wgcna_net"]]
str(wgcna_net)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 6, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
# 假设 adjMat 存在于 wgcna_net 的某个元素中，例如 "adjacency"
if ("adjacency" %in% names(wgcna_net)) {
  adjMat <- wgcna_net$adjacency
} else {
  # 如果不在 "adjacency"，检查其他可能的名称
  adjMat <- wgcna_net[["你的实际邻接矩阵名称"]]
}
datExpr <- seurat_obj@misc$datExpr
adjMat <- adjacency(datExpr, power = softPower)
# 确认 adjMat 的结构
print(dim(adjMat))
for (mod in 1:nrow(moduleColors)){ #【对每个模块进行循环】
  modules = names(table(moduleColors))[mod]  
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modProbes = probes[inModule]
  modGenes = modProbes  
  modTOM = TOM[inModule, inModule] 
  dimnames(modTOM) = list(modProbes, modProbes)
  edges_File = paste("CytoscapeInput-edges-", modules , ".txt", sep="")
  nodes_File = paste("CytoscapeInput-nodes-", modules, ".txt", sep="")
  outEdge=paste(cytoDir,edges_File,sep="\\")
  outNode=paste(cytoDir,nodes_File,sep="\\")
  cyt = exportNetworkToCytoscape(consTomDS,
                                 edgeFile = outEdge,
                                 nodeFile = outNode,
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
}
?exportNetworkToCytoscape
head(net)[1,]
### 不net### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)
plot_list <-ModuleFeaturePlot(
  seurat_obj,
  reduction = "umap",
  features = "hMEs",
  order_points = TRUE, # order so the points with highest hMEs are on top
  restrict_range = TRUE,
  point_size = 0.5,
  alpha = 1,
  label_legend = FALSE,
  raster_dpi = 500,
  raster_scale = 1,
  plot_ratio = 1,
  title = TRUE
)
plot_list[1]
wrap_plots(plot_list, ncol=4)#可以展示四种分数（hMEs, MEs, scores, or average）
ggsave(filename = '~/NSCLC/zong/fig/MALIGNANT/WGCNA_LUADplotUMAP.pdf',height = 5,width = 9,dpi = 300)
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)
str(wrap_plots(plot_list, ncol=6))  
ggsave(filename = '~/NSCLC/zong/fig/MALIGNANT/WGCNA_LUADplotUMAP-scores.pdf',height = 5,width = 9,dpi = 300)
ModuleCorrelogram(seurat_obj,
                  exclude_grey = TRUE, # 默认删除灰色模块
                  features = "hMEs") 
ggsave(filename = '~/NSCLC/zong/fig/MALIGNANT/WGCNA_LUADModuleCorrelogram.pdf',height = 5,width = 5,dpi = 300)
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
p <- DotPlot(seurat_obj, features=mods, group.by = 'Tissue')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
p
# plot output
ggsave( filename = '~/NSCLC/zong/fig/MALIGNANT/WGCNA_LUADDotplot-Tisue.pdf',height = 6,width = 6,dpi = 300,plot = p)
# 每个模块在不同样本中的情况
seurat_obj$leiden <- do.call(rbind, strsplit(as.character(seurat_obj$leiden), ' '))[,1]

p <- ModuleRadarPlot(
  seurat_obj,
  group.by = 'leiden',
  axis.label.size=4,
  grid.label.size=4
)
p
ggsave(filename = '~/NSCLC/zong/fig/MALIGNANT/WGCNA_LUADModuleRadarPlot-leiden.pdf',height = 8,width = 8,dpi = 300,plot = p)
options(future.globals.maxSize = 5 * 1024^3)  # 4 GB
dev.off()
pdf('~/NSCLC/zong/fig/MALIGNANT/WGCNA_LUSCHubGeneNetworkPlot.pdf',height = 8,width = 8)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
ggsave(filename = '~/NSCLC/zong/fig/MALIGNANT/WGCNA_LUADHubGeneNetworkPlot.pdf',height = 8,width = 8,dpi = 300)

# 定义enrichr databases
dbs <- c('GO_Biological_Process_2021',
         'GO_Cellular_Component_2021',
         'GO_Molecular_Function_2021')
# 富集分析
?RunEnrichr
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=c("GO_Biological_Process_2021", "GO_Cellular_Component_2021",
        "GO_Molecular_Function_2021"),
  max_genes = Inf# use max_genes = Inf to choose all genes
)
# 检索输出表
enrich_df <- GetEnrichrTable(seurat_obj)

# 查看结果
head(enrich_df)

# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more are shown if there are ties)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE# do you want to show the enrichment as a log scale?
)

# enrichr dotplot
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (default)
  database = "GO_Biological_Process_2023", # this must match one of the dbs used previously
  n_terms=2, # number of terms per module
  term_size=8, # font size for the terms
  p_adj = FALSE# show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))
