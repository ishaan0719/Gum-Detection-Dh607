library(Seurat)
library(Matrix)
library(edgeR)

seurat_obj <- readRDS(
  file.path("data", "processed", "GSE171213_HC_PD_annotated.rds")
)

# Use raw counts and sum them within each donor and cell type.
counts <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
metadata <- seurat_obj@meta.data
metadata$cell <- rownames(metadata)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

all_results <- list()
cell_types <- sort(unique(metadata$cell_type))

for (cell_type_name in cell_types) {
  cat("\nAnalyzing:", cell_type_name, "\n")

  cells_in_type <- metadata$cell[metadata$cell_type == cell_type_name]
  type_metadata <- metadata[cells_in_type, , drop = FALSE]

  # Sum cell counts for each donor.
  cells_by_donor <- split(rownames(type_metadata), type_metadata$donor)
  donor_counts <- vapply(
    cells_by_donor,
    function(cell_ids) {
      Matrix::rowSums(counts[, cell_ids, drop = FALSE])
    },
    FUN.VALUE = numeric(nrow(counts))
  )

  # Find the condition for each donor.
  sample_info <- unique(type_metadata[, c("donor", "condition")])
  sample_info <- sample_info[match(colnames(donor_counts), sample_info$donor), ]

  healthy_n <- sum(sample_info$condition == "Healthy")
  disease_n <- sum(sample_info$condition == "Periodontitis")

  if (healthy_n < 2 || disease_n < 2) {
    cat("Skipping: fewer than 2 donors in one condition.\n")
    next
  }

  group <- factor(
    sample_info$condition,
    levels = c("Healthy", "Periodontitis")
  )

  y <- DGEList(counts = donor_counts, group = group)
  keep_genes <- filterByExpr(y, group = group)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  design <- model.matrix(~ group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = 2)

  result <- topTags(test, n = Inf)$table
  result$gene <- rownames(result)
  result$cell_type <- cell_type_name
  result$healthy_donors <- healthy_n
  result$periodontitis_donors <- disease_n

  safe_name <- gsub("[^A-Za-z0-9]+", "_", cell_type_name)
  write.csv(
    result,
    file.path("results", "tables", paste0("DE_", safe_name, ".csv")),
    row.names = FALSE
  )

  all_results[[cell_type_name]] <- result
  cat("Tested genes:", nrow(result), "\n")
}

if (length(all_results) > 0) {
  combined <- do.call(rbind, all_results)
  write.csv(
    combined,
    file.path("results", "tables", "DE_all_cell_types.csv"),
    row.names = FALSE
  )
  cat("\nSaved combined results: results/tables/DE_all_cell_types.csv\n")
} else {
  cat("\nNo cell types had enough donors in both conditions to test.\n")
}