library(Seurat)
library(Matrix)
library(edgeR)
library(ggplot2)

obj <- readRDS(
  file.path("data", "processed", "GSE171213_HC_PD_annotated.rds")
)
counts <- LayerData(obj, assay = "RNA", layer = "counts")
meta <- obj@meta.data

candidates <- list(
  "Neutrophils" = c("LTF", "PGLYRP1", "CAMP"),
  "Epithelial cells" = c("HMOX1", "CCL20", "DST", "COL17A1", "CXCL6"),
  "Endothelial cells" = c("SELL", "THBS1", "CYR61", "MTHFR", "CRNDE")
)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

for (cell_type in names(candidates)) {
  type_meta <- meta[meta$cell_type == cell_type, , drop = FALSE]
  donors <- sort(unique(type_meta$donor))

  # Sum raw counts within each donor.
  donor_counts <- sapply(donors, function(donor) {
    cell_ids <- rownames(type_meta)[type_meta$donor == donor]
    Matrix::rowSums(counts[, cell_ids, drop = FALSE])
  })

  rownames(donor_counts) <- rownames(counts)

  # Convert donor-level counts to normalized log-CPM for plotting.
  dge <- DGEList(counts = donor_counts)
  dge <- calcNormFactors(dge)
  log_cpm <- cpm(dge, log = TRUE, prior.count = 2)

  genes_to_plot <- intersect(candidates[[cell_type]], rownames(log_cpm))
  if (length(genes_to_plot) == 0) {
    message("No candidate genes found for ", cell_type)
    next
  }

  plot_data <- expand.grid(
    gene = genes_to_plot,
    donor = colnames(log_cpm),
    stringsAsFactors = FALSE
  )
  plot_data$logCPM <- mapply(
    function(gene, donor) log_cpm[gene, donor],
    plot_data$gene,
    plot_data$donor
  )

  donor_conditions <- unique(type_meta[, c("donor", "condition")])
  plot_data <- merge(plot_data, donor_conditions, by = "donor")

  p <- ggplot(plot_data, aes(x = condition, y = logCPM, color = condition)) +
    geom_point(
      position = position_jitter(width = 0.08, height = 0),
      size = 3
    ) +
    facet_wrap(~ gene, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(
      title = paste("Donor-level expression:", cell_type),
      x = NULL,
      y = "Normalized log2 counts per million"
    ) +
    theme(legend.position = "none")

  filename <- gsub("[^A-Za-z0-9]+", "_", cell_type)
  ggsave(
    file.path("results", "figures", paste0("candidate_genes_", filename, ".png")),
    p,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}

cat("Saved candidate-gene donor plots in results/figures/.\n")