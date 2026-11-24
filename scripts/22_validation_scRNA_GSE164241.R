# Independent single-cell validation: GSE164241 (Williams et al., Cell 2021),
# healthy gingiva (GM, 13 samples) vs periodontitis gingiva (PD, 8 samples / 5 donors).
# QC -> clustering -> marker-score annotation into broad cell types -> donor-level
# composition, cell-type pseudobulk of the GSE171213 candidate genes, and
# tissue-level pseudobulk DE for genome-wide concordance.
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(edgeR)
  library(dplyr); library(ggplot2)})
set.seed(1)
indir <- "data/external/GSE164241"
mats <- list.files(indir, pattern = "_matrix.mtx.gz$")
objs <- list()
for (m in mats) {
  id <- sub("_matrix.mtx.gz$", "", m); sample <- sub("^GSM[0-9]+_", "", id)
  feat <- file.path(indir, paste0(id, "_features.tsv.gz"))
  if (!file.exists(feat)) feat <- file.path(indir, paste0(id, "_genes.tsv.gz"))
  x <- ReadMtx(mtx = file.path(indir, m), cells = file.path(indir, paste0(id, "_barcodes.tsv.gz")),
               features = feat, feature.column = 2)
  o <- CreateSeuratObject(x, project = sample, min.features = 200)
  o$sample <- sample
  o$donor <- sub("[a-z]$", "", sample)
  o$condition <- ifelse(grepl("^PD", sample), "Periodontitis", "Healthy")
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  o <- subset(o, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
  if (ncol(o) > 4000) o <- o[, sample(colnames(o), 4000)]   # random down-sampling keeps proportions
  cat(sample, ncol(o), "cells\n")
  objs[[sample]] <- o
}
obj <- merge(objs[[1]], objs[-1])
obj <- JoinLayers(obj)
obj <- NormalizeData(obj) |> FindVariableFeatures(nfeatures = 2000) |> ScaleData() |> RunPCA(npcs = 30)
obj <- FindNeighbors(obj, dims = 1:20) |> FindClusters(resolution = 0.8) |> RunUMAP(dims = 1:20)

markers <- list(
  "T/NK cells" = c("CD3D", "CD3E", "CD2", "NKG7", "GNLY"),
  "B cells" = c("MS4A1", "CD79A", "CD19", "BANK1"),
  "Plasma cells" = c("MZB1", "JCHAIN", "DERL3", "TNFRSF17"),
  "Neutrophils" = c("CSF3R", "FCGR3B", "CXCR2", "S100A8"),
  "Monocytes/macrophages" = c("CD68", "LYZ", "CD14", "C1QA", "AIF1"),
  "Mast cells" = c("TPSAB1", "CPA3", "KIT", "MS4A2"),
  "Endothelial cells" = c("PECAM1", "VWF", "CDH5", "ACKR1"),
  "Fibroblasts" = c("COL1A1", "DCN", "LUM", "COL3A1"),
  "Pericytes/smooth muscle" = c("RGS5", "ACTA2", "MYH11", "NOTCH3"),
  "Epithelial cells" = c("KRT14", "KRT5", "KRT13", "SFN"))
markers <- lapply(markers, intersect, rownames(obj))
obj <- AddModuleScore(obj, features = markers, name = "ms")
ms <- obj@meta.data[, paste0("ms", seq_along(markers))]; colnames(ms) <- names(markers)
cl_score <- aggregate(ms, by = list(cluster = obj$seurat_clusters), FUN = mean)
cl_label <- setNames(names(markers)[apply(cl_score[, -1], 1, which.max)], cl_score$cluster)
obj$cell_type <- unname(cl_label[as.character(obj$seurat_clusters)])
write.csv(cbind(cl_score, label = cl_label), "results/tables/GSE164241_cluster_annotation.csv", row.names = FALSE)

meta <- obj@meta.data; um <- Embeddings(obj, "umap")
write.csv(data.frame(cell = rownames(meta), sample = meta$sample, donor = meta$donor,
                     condition = meta$condition, cell_type = meta$cell_type,
                     UMAP_1 = um[, 1], UMAP_2 = um[, 2]),
          "results/tables/GSE164241_cells.csv", row.names = FALSE)

# Donor-level composition (multiple biopsies of one donor are pooled)
comp <- meta |> count(donor, condition, cell_type) |> group_by(donor) |>
  mutate(prop = n / sum(n)) |> ungroup()
comp <- tidyr::complete(comp, tidyr::nesting(donor, condition), cell_type,
                        fill = list(n = 0, prop = 0))
write.csv(comp, "results/tables/GSE164241_composition_by_donor.csv", row.names = FALSE)
tests <- comp |> group_by(cell_type) |> summarise(
  healthy_median = median(prop[condition == "Healthy"]) * 100,
  periodontitis_median = median(prop[condition == "Periodontitis"]) * 100,
  p_value = wilcox.test(prop[condition == "Healthy"], prop[condition == "Periodontitis"], exact = FALSE)$p.value) |>
  mutate(FDR_BH = p.adjust(p_value, "BH")) |> arrange(p_value)
write.csv(tests, "results/tables/GSE164241_composition_tests.csv", row.names = FALSE)
print(tests)

# Pseudobulk (donor x cell type) for candidate genes and tissue-level DE
counts <- LayerData(obj, assay = "RNA", layer = "counts")
run_pb <- function(cells_keep) {
  d <- meta$donor[cells_keep]
  pb <- sapply(split(which(cells_keep), d), function(ix) Matrix::rowSums(counts[, ix, drop = FALSE]))
  grp <- factor(ifelse(grepl("^PD", colnames(pb)), "Periodontitis", "Healthy"), c("Healthy", "Periodontitis"))
  if (min(table(grp)) < 2) return(NULL)
  y <- DGEList(pb, group = grp); y <- y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE]
  y <- normLibSizes(y); des <- model.matrix(~grp); y <- estimateDisp(y, des)
  r <- topTags(glmQLFTest(glmQLFit(y, des, robust = TRUE), coef = 2), n = Inf)$table
  r$gene <- rownames(r); r
}
tissue <- run_pb(rep(TRUE, ncol(obj)))
write.csv(tissue, "results/tables/GSE164241_DE_tissue_pseudobulk.csv", row.names = FALSE)
res <- list()
for (ct in sort(unique(meta$cell_type))) {
  r <- run_pb(meta$cell_type == ct); if (is.null(r)) next
  r$cell_type <- ct; res[[ct]] <- r
}
write.csv(do.call(rbind, res), "results/tables/GSE164241_DE_celltype_pseudobulk.csv", row.names = FALSE)
cat("Cells:", ncol(obj), " tissue DE FDR<0.05:", sum(tissue$FDR < 0.05), "\n")

p <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE) + ggtitle("GSE164241: broad cell types") + theme_bw()
ggsave("results/figures/GSE164241_umap_cell_types.png", p, width = 9, height = 7, dpi = 200, bg = "white")
p2 <- DotPlot(obj, features = unique(unlist(lapply(markers, head, 3))), group.by = "cell_type") + RotatedAxis()
ggsave("results/figures/GSE164241_marker_dotplot.png", p2, width = 13, height = 5.5, dpi = 200, bg = "white")
