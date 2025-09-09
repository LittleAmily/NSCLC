library(Seurat)
library(hdWGCNA)
library(ggplot2)
library(patchwork)
setwd("/home/data/sdzl14")
LUSC <- readRDS(file = '~/NSCLC/zong/LUSC.rds')
LUAD  <- readRDS(file = '~/NSCLC/zong/LUAD_qc.rds')
merged_object <- merge(x = LUSC, y = LUAD,  add.cell.ids = c("LUSC", "LUAD"))
Idents(merged_object) <- 'Tissue'

source('/home/data/sdzl14/NSCLC/code/run_hdWGCNA_for_tissue.r')

# 获取所有 tissue 类型
tissues <- unique(merged_object$Tissue)

# 设置输出目录
output_base_dir <- "/NSCLC/zong/fig/MALIGNANT/WCGNA1"

# 创建主输出文件夹（如果不存在）
dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

# 对每个 tissue 执行流程
for (tissue in tissues) {
  tissue_output_dir <- file.path(output_base_dir, tissue)
  dir.create(tissue_output_dir, showWarnings = FALSE, recursive = TRUE)
  run_hdWGCNA_for_tissue(tissue, merged_object, tissue_output_dir)
}

tissue_name <- "tumor_metastasis"

print(paste0("Processing tissue: ", tissue_name))

subset <- subset(merged_object, ident = tissue_name)
subset <- NormalizeData(subset)
subset <- ScaleData(subset)
subset <- FindVariableFeatures(subset)
subset <- RunPCA(subset)
subset <- RunUMAP(subset, dims = 1:20)
seurat_obj <- SetupForWGCNA(
  subset,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = tissue_name
)
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Timepoint", "Tissue", "Origin", "Pathtype", "leiden"),
  reduction = 'umap',
  k = 25,
  max_shared = 10,
  ident.group = 'Tissue'
)
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = 'Pre',
  group.by = 'Timepoint',
  assay = 'RNA',
  slot = 'counts'
)
seurat_obj <- TestSoftPowers(seurat_obj, networkType = 'signed')

wrap_plots(plot_list)
ggsave(filename = paste0("~/NSCLC/zong/fig/MALIGNANT/WGCNA_", tissue_name, "_PlotSoftPowers.pdf"), 
        dpi = 300, height = 10, width = 10)
seurat_obj <- ConstructNetwork(
  seurat_obj,
  setDatExpr = FALSE,
  corType = "pearson",
  networkType = "signed",
  TOMType = "signed",
  detectCutHeight = 0.995,
  minModuleSize = 25,
  mergeCutHeight = 0.1,
  overwrite_tom = TRUE,
  tom_name = tissue_name
)

seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(seurat_obj, scale.model.use = "linear", pc_dim = 1)
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'Timepoint', group_name = 'Pre',
  corFnc = "bicor", corOptions = "use='p'",
  harmonized = TRUE
)
output_dir <- "~/NSCLC/zong/fig/MALIGNANT"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(seurat_obj, file = file.path(output_dir, paste0(tissue_name, "_hdWGCNA_object.rds")))
dev.off()
pdf(file.path(output_dir, paste0("WGCNA_", tissue_name, "_PlotDendrogram.pdf")), 
    height = 8, width = 8)

PlotDendrogram(seurat_obj, main = paste0(tissue_name, " hdWGCNA Dendrogram"))
dev.off()
