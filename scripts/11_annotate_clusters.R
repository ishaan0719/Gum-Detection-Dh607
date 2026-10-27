library(Seurat)
library(ggplot2)

seurat_obj <- readRDS(
  file.path("data", "processed", "GSE171213_HC_PD_seurat.rds")
)

cluster_labels <- c(
  "0"  = "T cells",
  "1"  = "Cytotoxic T cells",
  "2"  = "NK/cytotoxic lymphocytes",
  "3"  = "Endothelial cells",
  "4"  = "Plasma cells",
  "5"  = "Neutrophils",
  "6"  = "B cells",
  "7"  = "Monocytes/macrophages",
  "8"  = "T cells",
  "9"  = "Neutrophils",
  "10" = "Fibroblasts",
  "11" = "Mast cells",
  "12" = "Plasma cells",
  "13" = "Epithelial cells",
  "14" = "Proliferating cells",
  "15" = "Pericytes/smooth-muscle-like",
  "16" = "Uncertain myeloid-like",
  "17" = "Neutrophils",
  "18" = "Erythroid cells"
)

cluster_ids <- as.character(seurat_obj$seurat_clusters)
seurat_obj$cell_type <- unname(cluster_labels[cluster_ids])

if (any(is.na(seurat_obj$cell_type))) {
  stop("At least one cluster is missing a label.")
}

p <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) +
  theme_bw() +
  ggtitle("Provisional cell-type labels")

ggsave(
  "results/figures/umap_cell_types.png",
  p,
  width = 11,
  height = 8,
  dpi = 300,
  bg = "white"
)

saveRDS(
  seurat_obj,
  file.path("data", "processed", "GSE171213_HC_PD_annotated.rds")
)

cat("Saved provisional labels and UMAP plot.\n")