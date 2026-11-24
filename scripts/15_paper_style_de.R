library(Seurat)

obj <- readRDS(
  "data/processed/GSE171213_HC_PD_annotated.rds"
)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

all_results <- list()
cell_types <- sort(unique(obj$cell_type))

for (ct in cell_types) {
  message("Analyzing: ", ct)

  cells <- colnames(obj)[obj$cell_type == ct]
  sub_obj <- subset(obj, cells = cells)

  # Check that this cell type is represented by at least two donors
  # in both conditions.
  donor_info <- unique(
    data.frame(
      donor = sub_obj$donor,
      condition = sub_obj$condition
    )
  )

  donors_per_condition <- table(donor_info$condition)

  if (!all(c("Healthy", "Periodontitis") %in% names(donors_per_condition)) ||
      donors_per_condition["Healthy"] < 2 ||
      donors_per_condition["Periodontitis"] < 2) {
    message("Skipping: too few donors in one condition.")
    next
  }

  Idents(sub_obj) <- "condition"

  result <- FindMarkers(
    sub_obj,
    ident.1 = "Periodontitis",
    ident.2 = "Healthy",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.10,
    only.pos = FALSE,
    return.thresh = 1
  )

  result$gene <- rownames(result)
  result$cell_type <- ct

  # Add a BH-adjusted p-value for comparison.
  result$FDR_BH <- p.adjust(result$p_val, method = "BH")

  all_results[[ct]] <- result
}

if (length(all_results) == 0) {
  stop("No cell types had enough donors in both conditions.")
}

combined <- do.call(rbind, all_results)
write.csv(
  combined,
  "results/tables/DE_paper_style_wilcoxon.csv",
  row.names = FALSE
)

message("Saved results/tables/DE_paper_style_wilcoxon.csv")