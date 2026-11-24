# Export log-normalised expression of variable genes, UMAP coordinates and labels
# from the annotated GSE171213 object for the Python ML / plotting steps, plus
# donor x cell-type pseudobulk counts.
library(Seurat)
library(Matrix)

obj <- readRDS("data/processed/GSE171213_HC_PD_annotated.rds")
out <- "data/processed/ml_input"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

hvg <- VariableFeatures(obj)
data <- LayerData(obj, assay = "RNA", layer = "data")[hvg, ]
writeMM(t(data), file.path(out, "lognorm_hvg_cells_x_genes.mtx"))
writeLines(hvg, file.path(out, "hvg.txt"))

meta <- obj@meta.data
umap <- Embeddings(obj, "umap")
meta$UMAP_1 <- umap[, 1]; meta$UMAP_2 <- umap[, 2]
meta$cell <- rownames(meta)
write.csv(meta, file.path(out, "cell_metadata.csv"), row.names = FALSE)

# Donor-level pseudobulk across all cells (for donor-level classification).
counts <- LayerData(obj, assay = "RNA", layer = "counts")
donors <- sort(unique(meta$donor))
pb <- sapply(donors, function(d) Matrix::rowSums(counts[, meta$donor == d, drop = FALSE]))
write.csv(pb, file.path(out, "pseudobulk_donor_counts.csv"))
cat("Exported", ncol(data), "cells x", length(hvg), "genes\n")
