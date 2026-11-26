# Project-GumDetective

Single-cell and bulk transcriptomics of healthy gingiva vs periodontitis.

## Data
- `data/raw/GSE171213_AllSample.counts.tsv.gz`: discovery scRNA-seq (GEO GSE171213)
- `data/external/`: downloaded by the scripts below. GSE10334 and GSE16134 series matrices, GPL570 annotation, GSE164241 GM/PD 10x files (from GEO).

## Pipeline (run from the project root, in order)
| Step | Script | What it does |
|---|---|---|
| 01–08 | `*.R / *.py` | inspect, convert to sparse, select HC1–4/PD1–5, QC metrics, filter, export |
| 09–11 | Seurat | normalise, cluster, UMAP, markers, annotation |
| 12 | `12_pseudobulk_de.R` | donor × cell-type edgeR pseudobulk DE |
| 13 | candidate gene plots | |
| 14–15 | paper-style QC and Wilcoxon DE | reproduction of Chen et al. |
| 16 | composition per donor | |
| 17 | `17_prepare_bulk_validation.py` | bulk cohorts, probe→gene, GSE16134 duplicate detection |
| 18 | `18_export_for_ml.R` | export HVG matrix / metadata / pseudobulk |
| 19 | `19_pathway_enrichment.py` | Enrichr GO/KEGG per cell type |
| 20 | `20_tissue_pseudobulk_de.R` | tissue-level pseudobulk DE |
| 21 | `21_bulk_validation_and_ml.py` | bulk DE, signature tests, hypothesis-testing vs LR/SVM/RF, external test |
| 22 | `22_validation_scRNA_GSE164241.R` | independent scRNA cohort |
| 23 | `23_cross_dataset_and_scRNA_ml.py` | cross-cohort replication, single-cell leakage analysis |
| 24 | `24_make_figures.py` | final figures → `results/figures/final` |

Environment: R 4.x with Seurat 5, edgeR; Python 3 with pandas, scipy, scikit-learn, statsmodels, matplotlib, gseapy.
