<div align="center">

# 🦷 GumDetective

### Periodontitis at single-cell resolution: what changes, and what replicates?

**Ishaan Srivastava** (23B2187) · **Neha Kanwadiya** (23B1852)

![R](https://img.shields.io/badge/R-4.x-276DC3?logo=r&logoColor=white)
![Seurat](https://img.shields.io/badge/Seurat-5-1C7C78)
![Python](https://img.shields.io/badge/Python-3-3776AB?logo=python&logoColor=white)
![scikit-learn](https://img.shields.io/badge/scikit--learn-ML-F7931E?logo=scikitlearn&logoColor=white)
![Data](https://img.shields.io/badge/data-GEO%20(4%20cohorts)-B5503F)

</div>

<p align="center">
  <img src="docs/images/normal_vs_periodontitis.webp" alt="Healthy gums and tooth compared with periodontitis: red, inflamed gums, gum recession, and plaque and tartar buildup" width="720">
</p>

**Periodontitis** is a chronic inflammatory gum disease that destroys the tissue and bone holding the teeth. We used single-cell and bulk RNA data to find which cells and genes change in diseased gums, then tested every finding in three independent cohorts. Throughout, the **donor or patient is the unit of analysis, never the individual cell**.

---

## 🎯 Hypothesis and verdict

> *Periodontitis is associated with changes in periodontal cell populations and gene expression.*

| Part | Verdict | Evidence |
|---|:---:|---|
| Gene expression changes | ✅ **Supported** | 9,929 genes change in 247 biopsies · 39 genes up in all 3 cohorts · AUC 0.93 |
| Cell populations change | 🟡 **Partly** | Plasma cells ↑ and macrophages ↑ in both cohorts; four other shifts reverse |
| scRNA signature works in bulk | 🔴 **Weak** | AUC 0.63, not confirmed externally |

## 🔬 Key results

| # | Finding |
|:---:|---|
| 1 | **Cells are not replicates:** cell-level testing gives 4,054 hits; donor-level pseudobulk gives **14** |
| 2 | **Cells:** plasma cells and macrophages rise in both single-cell cohorts; endothelial and T-cell trends reverse |
| 3 | **Genes:** **THBS1, CYR61, CXCL6** replicate in every cohort; neutrophil genes do not |
| 4 | **Pathways:** inflammation up, epithelial structure down; NF-κB and osteoclast pathways in neutrophils |
| 5 | **Bulk tissue:** strong signal (9,929 genes), but the scRNA signature alone is weak (AUC 0.63) |
| 6 | **Tests vs ML:** simple DE score **0.91** vs logistic regression/SVM **0.93**; 0.90–0.99 on 29 new patients |
| 7 | **Leakage:** single-cell classifier drops from **0.94** (random split) to **0.50** on unseen donors |

## 📦 Data

All public on GEO; not stored here (`data/` is git-ignored).

| Role | Dataset | Type | Size |
|---|---|---|---|
| Discovery | [GSE171213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171213) | scRNA-seq | 4 healthy vs 5 periodontitis donors · 27,272 cells |
| Validation | [GSE164241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164241) | scRNA-seq | 13 healthy vs 5 periodontitis donors |
| Validation + ML | [GSE10334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10334) | Microarray | 247 biopsies · 90 patients |
| External test | [GSE16134](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16134) | Microarray | 63 samples · 29 new patients |

> ⚠️ GSE16134 re-uploads **245** GSE10334 arrays under new IDs. Script 17 detects and removes them so the external test is truly independent.

<details>
<summary>Where to put the downloaded files</summary>

```
data/raw/GSE171213_AllSample.counts.tsv.gz
data/external/GSE10334_series_matrix.txt.gz
data/external/GSE16134_series_matrix.txt.gz
data/external/GPL570.annot.gz
data/external/GSE164241/   # 10x files: *_matrix.mtx.gz, *_barcodes.tsv.gz, *_features.tsv.gz
```
</details>

## ⚙️ Pipeline

Run `scripts/01` → `24` in order from the project root.

```
QC & filtering  →  Clustering & labels  →  Donor-level DE  →  Replication  →  Pathways  →  Tests, then ML
   (01–08)            (09–11)              (12–16, 20)       (17, 22–23)      (19)          (21, 23)
```

| Scripts | Step |
|---|---|
| 01–08 | Select HC1–4 / PD1–5, QC (>200 genes, <25% mito), filter, export |
| 09–11 | Seurat: normalise, PCA, UMAP, clustering, cell-type labels |
| 12, 20 | edgeR pseudobulk DE per cell type and per tissue |
| 13–16 | Candidate plots, paper-style Wilcoxon reproduction, composition tests |
| 17, 22 | Bulk cohort parsing + duplicate check; GSE164241 validation |
| 18–19 | ML export; Enrichr GO/KEGG enrichment |
| 21, 23 | Bulk tests vs LR/SVM/RF (patient-grouped CV), cross-cohort replication, leakage test |
| 24 | Final figures |

## 📊 Results in this repository

Everything is committed, so you can read the findings without running anything.

**Figures** · `results/figures/final/`

| Figure | Shows |
|---|---|
| `fig1_umap` | Cell types and condition on UMAP |
| `fig2_composition_two_cohorts` | Per-donor cell proportions, both cohorts |
| `fig3_de_method_comparison` | Cell-level vs donor-level DE hits |
| `fig4_candidate_replication` | 14 candidate genes across 3 cohorts |
| `fig5_enrichment` | Enriched pathways by cell type |
| `fig6_bulk_volcano` | Bulk biopsies volcano plot |
| `fig7_ml_bulk_comparison` | Tests vs ML performance |
| `fig8_roc_bulk` | ROC curves |
| `fig9_leakage` | Random split vs held-out donors |

**Key tables** · `results/tables/`

| Table | Contents |
|---|---|
| `DE_all_cell_types.csv` | Donor-level DE, all cell types |
| `composition_two_scRNA_cohorts.csv` | Cell composition tests |
| `candidate_genes_replication.csv` | Candidate genes in 3 cohorts |
| `replicated_genes_three_cohorts.csv` | The 39 shared genes |
| `signature_validation_bulk.csv` | scRNA signatures tested in bulk |
| `ml_bulk_cv_summary.csv` · `ml_bulk_external_test.csv` | ML results |

## 🧰 Environment, limitations and next steps

**Environment:** R 4.x (Seurat 5, edgeR) · Python 3 (pandas, scipy, scikit-learn, statsmodels, matplotlib, gseapy). Single-cell steps need plenty of memory; Enrichr needs internet.

**Limitations:** small discovery cohort (4 vs 5 donors) · only 5 healthy samples in the external test · cohorts differ in sampling and technology · no dedicated ambient-RNA/doublet removal.

**Next steps:** ambient-RNA and doublet removal · more donors · analyse post-treatment donors · protein-level validation of THBS1, CYR61, CXCL6.
