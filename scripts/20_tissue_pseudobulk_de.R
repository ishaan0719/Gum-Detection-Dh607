# Tissue-level (all cells per donor) pseudobulk DE, periodontitis vs healthy, edgeR QL.
# This is the single-cell analogue of a bulk biopsy and is used to derive a
# signature that is then tested in the independent bulk cohorts.
library(edgeR)
pb <- as.matrix(read.csv("data/processed/ml_input/pseudobulk_donor_counts.csv",
                         row.names = 1, check.names = FALSE))
group <- factor(ifelse(grepl("^HC", colnames(pb)), "Healthy", "Periodontitis"),
                levels = c("Healthy", "Periodontitis"))
y <- DGEList(pb, group = group)
y <- y[filterByExpr(y, group = group), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
res <- topTags(glmQLFTest(fit, coef = 2), n = Inf)$table
res$gene <- rownames(res)
write.csv(res, "results/tables/DE_tissue_pseudobulk.csv", row.names = FALSE)
write.csv(cpm(y, log = TRUE, prior.count = 2), "data/processed/ml_input/pseudobulk_donor_logcpm.csv")
cat("Genes tested:", nrow(res), " FDR<0.05:", sum(res$FDR < 0.05),
    " FDR<0.1:", sum(res$FDR < 0.1), " p<0.001:", sum(res$PValue < 0.001), "\n")
print(head(res[, c("gene", "logFC", "PValue", "FDR")], 25))
