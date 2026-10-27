library(Matrix)
library(Seurat)
library(ggplot2)

input_dir <- file.path("data", "processed", "seurat_input")
results_dir <- file.path("results", "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Read the sparse counts matrix (genes × cells).
counts <- readMM(gzfile(file.path(input_dir, "counts.mtx.gz")))
genes <- readLines(file.path(input_dir, "genes.txt"))
cells <- readLines(file.path(input_dir, "cells.txt"))

if (nrow(counts) != length(genes)) {
  stop("Gene names do not match the matrix rows.")
}
if (ncol(counts) != length(cells)) {
  stop("Cell names do not match the matrix columns.")
}

rownames(counts) <- make.unique(genes)
colnames(counts) <- cells

# Attach donor and condition labels.
metadata <- read.csv(
  file.path(input_dir, "metadata.csv"),
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$cell
metadata$cell <- NULL
metadata <- metadata[colnames(counts), , drop = FALSE]

if (!identical(rownames(metadata), colnames(counts))) {
  stop("Metadata cell order does not match the count matrix.")
}

# Create Seurat object and calculate mitochondrial percentage.
seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  min.cells = 3
)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)

# Normalize, identify variable genes, and run PCA.
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# Cluster cells and calculate a UMAP embedding.
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Save the object for the next analysis steps.
saveRDS(
  seurat_obj,
  file = file.path("data", "processed", "GSE171213_HC_PD_seurat.rds")
)

# Save UMAP plots colored by cluster and condition.
p_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
  ggtitle("GSE171213: clusters")

p_condition <- DimPlot(seurat_obj, reduction = "umap", group.by = "condition") +
  ggtitle("GSE171213: condition")

ggsave(
  file.path(results_dir, "umap_clusters.png"),
  p_clusters,
  width = 8, height = 6, dpi = 300
)
ggsave(
  file.path(results_dir, "umap_condition.png"),
  p_condition,
  width = 8, height = 6, dpi = 300
)

cat("Cells:", ncol(seurat_obj), "\n")
cat("Genes:", nrow(seurat_obj), "\n")
cat("Clusters:", length(levels(Idents(seurat_obj))), "\n")
cat("Saved Seurat object and UMAP plots.\n")