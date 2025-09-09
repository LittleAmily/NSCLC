# 定义主函数
run_hdWGCNA_for_tissue <- function(tissue_name, data, output_dir) {
  print(paste0("Processing tissue: ", tissue_name))
  Idents(data) <- "Tissue"
  # Step 1: Subset and preprocess
  subset <- subset(data, ident = tissue_name)
  subset <- NormalizeData(subset)
  subset <- ScaleData(subset)
  subset <- FindVariableFeatures(subset)
  subset <- RunPCA(subset)
  subset <- RunUMAP(subset, dims = 1:20)
  # Step 2: Setup for WGCNA
  seurat_obj <- SetupForWGCNA(
    subset,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = tissue_name
  )
  
  # Step 3: Metacells by groups
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("Timepoint", "Tissue", "Origin", "Pathtype", "leiden"),
    reduction = 'umap',
    k = 25,
    max_shared = 10,
    ident.group = 'Tissue'
  )
  
  seurat_obj <- NormalizeMetacells(seurat_obj)
  
  # Step 4: SetDatExpr
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = 'Pre',
    group.by = 'Timepoint',
    assay = 'RNA',
    slot = 'counts'
  )
  
  # Step 5: Test soft powers
  seurat_obj <- TestSoftPowers(seurat_obj, networkType = 'signed')
  plot_list <- PlotSoftPowers(seurat_obj)
  wrap_plots(plot_list)
  ggsave(filename = file.path(output_dir, paste0("WGCNA_", tissue_name, "_PlotSoftPowers.pdf")),
          dpi = 300, height = 10, width = 10)
  
  # Step 6: Construct Network


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
  
  # Step 7: Plot Dendrogram
  # Step 8: Module Eigengenes & Connectivity
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
  seurat_obj <- ModuleEigengenes(seurat_obj, scale.model.use = "linear", pc_dim = 1)
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = 'Timepoint', group_name = 'Pre',
    corFnc = "bicor", corOptions = "use='p'",
    harmonized = TRUE
  )
  
  # Step 9: Save object
  saveRDS(seurat_obj, file = file.path(output_dir, paste0(tissue_name, "_hdWGCNA_object.rds")))
  
  # Step 10: Visualization
  
  dev.off()
  options(future.globals.maxSize = 5 * 1024^3)  # Increase limit if needed

  pdf(file.path(output_dir, paste0("WGCNA_", tissue_name, "_HubGeneNetworkPlot.pdf")), height = 8, width = 8)
  HubGeneNetworkPlot(seurat_obj, n_hubs = 3, n_other = 5, edge_prop = 0.75, mods = 'all')
  dev.off()
  
  print(paste0("Finished processing tissue: ", tissue_name))
}