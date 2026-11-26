# GumDetective

**Periodontitis at single-cell resolution: what changes, and what replicates?**

Ishaan Srivastava (23B2187) · Neha Kanwadiya (23B1852)

GumDetective compares healthy gum tissue (gingiva) with periodontitis using single-cell and bulk transcriptomics. We find disease-related changes in one single-cell dataset, then test whether they hold up in three independent cohorts. Every test treats the **donor or patient as the unit, never the individual cell**.

## Hypothesis

> Periodontitis is associated with changes in periodontal cell populations and gene expression.

**Verdict:** supported for gene expression, and partly supported for cell populations.

| Part of the hypothesis | Verdict | Evidence |
|---|---|---|
| Gene expression changes | Supported | 9,929 genes change in 247 biopsies; 39 genes up in all three cohorts; THBS1, CYR61, CXCL6 replicate; classifiers reach AUC 0.93 (0.90–0.99 on new patients) |
| Cell populations change | Partly supported | No significant shift in the 9-donor discovery cohort; plasma cells and macrophages rise in both single-cell cohorts (significant in the larger one); endothelial, pericyte, epithelial and T-cell shifts reverse |
| Single-cell signature carries over to bulk tissue | Weak | Significant but small (AUC 0.63); not confirmed in the external cohort |

## Key results

1. **Cells are not replicates.** The original paper's cell-level Wilcoxon test gives 4,054 significant gene–cell-type hits (2,389 unique genes). Donor-level pseudobulk with edgeR gives 14.
2. **Cell composition.** Plasma cells and monocytes/macrophages increase in both single-cell cohorts (FDR 0.02 in GSE164241). The endothelial increase and T-cell decrease seen in the discovery cohort reverse in the validation cohort.
3. **Candidate genes.** Of the 14 donor-level genes, THBS1, CYR61 and CXCL6 are up and significant in every cohort. SELL matches in direction only (p = 0.08 in GSE164241). The neutrophil genes LTF, CAMP and PGLYRP1 do not replicate.
4. **Pathways.** In bulk tissue, inflammatory and chemokine signalling go up and epidermis/keratinocyte genes go down. In the discovery cohort, neutrophils show NF-κB and osteoclast-differentiation pathways, and plasma cells show antibody-secretion machinery (single cohort, so interpretation only).
5. **Bulk tissue.** 9,929 genes differ between diseased and healthy sites in GSE10334 (FDR < 0.05); 39 genes are up in all three cohorts and none is consistently down. The single-cell tissue signature separates sites only weakly (AUC 0.63, p = 0.002; 0.61 in GSE16134).
6. **Tests before ML.** A 50-gene DE score reaches AUC 0.91 under patient-grouped cross-validation; logistic regression and linear SVM reach 0.93. On 29 new patients the models reach 0.90–0.99, with only 5 healthy samples in that set.
7. **Donor leakage.** A cell-level classifier scores AUC 0.94 with a random cell split but 0.50 when each donor is held out, so it was recognising donors, not disease.

## Data

All data are public on GEO. They are not stored in this repository (`data/` is git-ignored).

| Role | Dataset | Type | Samples |
|---|---|---|---|
| Discovery | [GSE171213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171213) | scRNA-seq | 4 healthy (HC1–4) vs 5 untreated periodontitis (PD1–5) donors; 27,272 cells after QC (from 45,191); 3 post-treatment donors excluded |
| Validation | [GSE164241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164241) | scRNA-seq (10x) | 13 healthy vs 5 periodontitis donors |
| Validation + ML training | [GSE10334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10334) | Microarray, GPL570 | 247 biopsies from 90 patients |
| External test | [GSE16134](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16134) | Microarray, GPL570 | 63 independent samples from 29 patients |

**Duplicate check.** GSE16134 re-uploads 245 GSE10334 arrays under new sample IDs. Script 17 finds them by correlating rank-transformed expression profiles and removes every affected patient, so the external test shares no samples with the training data.

Place the downloaded files as follows before running:

```
data/raw/GSE171213_AllSample.counts.tsv.gz
data/external/GSE10334_series_matrix.txt.gz
data/external/GSE16134_series_matrix.txt.gz
data/external/GPL570.annot.gz
data/external/GSE164241/   # the GM/PD 10x files: *_matrix.mtx.gz, *_barcodes.tsv.gz, *_features.tsv.gz (or *_genes.tsv.gz)
```

## Pipeline

Run the scripts from the project root, in order.

| Step | Script(s) | What it does |
|---|---|---|
| 01–02 | `01_inspect_data.R`, `02_preview_counts.R` | Inspect the GEO files and preview the count matrix |
| 03–04 | `03_convert_counts.py`, `04_select_primary_samples.py` | Convert counts to a sparse matrix; keep HC1–4 and PD1–5 |
| 05–08 | `05_calculate_qc.py` … `08_export_for_seurat.py` | QC metrics per cell and donor; filter (>200 genes, <25% mitochondrial, donor-specific upper limits); export for Seurat |
| 09–11 | `09_seurat_clustering.R` … `11_annotate_clusters.R` | Normalise, 2,000 variable genes, PCA (30), UMAP and clustering (20 PCs, resolution 0.5), marker genes, cell-type labels |
| 12 | `12_pseudobulk_de.R` | Donor × cell-type pseudobulk DE with edgeR (quasi-likelihood F-test) |
| 13 | `13_plot_candidate_genes.R` | Plots of the candidate genes |
| 14–15 | `14_compare_paper_qc.py`, `15_paper_style_de.R` | Reproduce the paper's QC and cell-level Wilcoxon DE (Chen et al.) |
| 16 | `16_cell_type_composition.R` | Cell-type proportions per donor; Mann–Whitney tests with BH FDR |
| 17 | `17_prepare_bulk_validation.py` | Parse bulk cohorts, collapse probes to genes, detect GSE16134 duplicates |
| 18 | `18_export_for_ml.R` | Export HVG matrix, metadata and pseudobulk for ML |
| 19 | `19_pathway_enrichment.py` | Enrichr GO/KEGG enrichment per cell type and in bulk |
| 20 | `20_tissue_pseudobulk_de.R` | Tissue-level pseudobulk DE |
| 21 | `21_bulk_validation_and_ml.py` | Bulk DE, signature tests, hypothesis tests vs logistic regression, SVM and random forest (patient-grouped 5-fold CV × 5), external test |
| 22 | `22_validation_scRNA_GSE164241.R` | Same QC, labelling, composition and pseudobulk DE on GSE164241 |
| 23 | `23_cross_dataset_and_scRNA_ml.py` | Cross-cohort replication; single-cell classifier with random split vs leave-one-donor-out |
| 24 | `24_make_figures.py` | Final figures in `results/figures/final` |

## Results in this repository

All outputs are committed, so the findings can be read without re-running the pipeline.

**Figures** (`results/figures/final`)

| File | Shows |
|---|---|
| `fig1_umap.png` | UMAP of annotated cell types and condition (GSE171213) |
| `fig2_composition_two_cohorts.png` | Per-donor cell-type proportions in both single-cell cohorts |
| `fig3_de_method_comparison.png` | Cell-level Wilcoxon vs donor-level pseudobulk hits per cell type |
| `fig4_candidate_replication.png` | Heatmap of the 14 candidate genes across three cohorts |
| `fig5_enrichment.png` | Enriched GO/KEGG terms by cell type |
| `fig6_bulk_volcano.png` | Volcano plot of GSE10334 biopsies |
| `fig7_ml_bulk_comparison.png` | ROC-AUC and balanced accuracy: tests vs ML, with external test |
| `fig8_roc_bulk.png` | Out-of-fold ROC curves in GSE10334 |
| `fig9_leakage.png` | Cell-level classifier under random vs donor-held-out validation |

**Main tables** (`results/tables`)

| File | Contents |
|---|---|
| `DE_all_cell_types.csv` | Donor-level pseudobulk DE for every cell type |
| `DE_paper_style_wilcoxon.csv` | Reproduced cell-level Wilcoxon DE |
| `cell_type_proportion_tests.csv`, `composition_two_scRNA_cohorts.csv` | Composition tests in one and both single-cell cohorts |
| `candidate_genes_replication.csv` | The 14 candidate genes in all three cohorts |
| `replicated_genes_three_cohorts.csv` | The 39 genes up in all three cohorts |
| `DE_bulk_GSE10334.csv`, `enrichment_bulk_GSE10334.csv` | Bulk DE and enrichment |
| `signature_validation_bulk.csv` | Single-cell gene sets tested in bulk |
| `ml_bulk_cv_summary.csv`, `ml_bulk_external_test.csv` | Cross-validation and external-test performance |
| `ml_scRNA_LODO_donor_predictions.csv` | Leave-one-donor-out predictions per donor |
| `bulk_summary.json`, `cross_dataset_summary.json` | Headline numbers |

## Environment

- R 4.x with Seurat 5 and edgeR
- Python 3 with pandas, scipy, scikit-learn, statsmodels, matplotlib and gseapy (Enrichr queries need internet access)

The single-cell steps need substantial memory. Script 22 down-samples GSE164241 to at most 4,000 cells per sample, which keeps cell-type proportions.

## Limitations

- The discovery cohort is small (4 vs 5 donors), so moderate composition shifts cannot be detected there.
- The external test set holds only 5 healthy samples, so its AUCs carry wide confidence intervals.
- The cohorts differ in sampling, technology and annotation depth.
- Ambient RNA and doublets are not removed beyond the QC thresholds.
- Single-cell pathway results come from one cohort.

**Next steps:** ambient-RNA and doublet removal, more donors, analysis of the three post-treatment donors, and protein-level validation of THBS1, CYR61 and CXCL6.
