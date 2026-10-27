library(Seurat)
library(ggplot2)

seurat_obj <- readRDS(
  file.path("data", "processed", "GSE171213_HC_PD_seurat.rds")
)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Find positive marker genes for each cluster.
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers,
  "results/tables/cluster_markers_all.csv",
  row.names = FALSE
)

# Save the top 10 markers per cluster for easier review.
top_markers <- do.call(
  rbind,
  lapply(
    split(markers, markers$cluster),
    function(x) head(x[order(x$p_val_adj), ], 10)
  )
)

write.csv(
  top_markers,
  "results/tables/cluster_markers_top10.csv",
  row.names = FALSE
)

# Plot selected reference markers to help annotate cell types.
reference_markers <- c(
  "PTPRC",                 # immune cells
  "CD3D", "CD3E",          # T cells
  "CD79A", "MS4A1",        # B cells
  "MZB1", "SDC1",          # plasma cells
  "FCGR3B", "CSF3R",       # neutrophils
  "LYZ", "S100A8",         # myeloid cells
  "COL1A1", "DCN",         # fibroblasts
  "PECAM1", "VWF",         # endothelial cells
  "EPCAM", "KRT19",         # epithelial cells
  "TPSAB1"                 # mast cells
)

reference_markers <- intersect(
  reference_markers,
  rownames(seurat_obj)
)

if (length(reference_markers) > 0) {
  p <- DotPlot(
    seurat_obj,
    features = reference_markers
  ) +
    RotatedAxis() +
    ggtitle("Reference marker expression by cluster")

  ggsave(
    "results/figures/reference_marker_dotplot.png",
    p,
    width = 12,
    height = 7,
    dpi = 300
  )
}

cat("Saved all cluster markers and top 10 markers per cluster.\n")
cat("Clusters:", length(unique(markers$cluster)), "\n")
cat("Reference marker plot saved if matching genes were present.\n")
