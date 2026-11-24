"""Final summary figures (results/figures/final)."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

OUT = "results/figures/final"; T = "results/tables"
HC, PD = "#2a78d6", "#eb6834"; INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
plt.rcParams.update({"font.family": "Helvetica", "font.size": 11, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "figure.facecolor": "white", "savefig.dpi": 220, "savefig.bbox": "tight"})
def save(fig, name): fig.savefig(f"{OUT}/{name}.png"); plt.close(fig)

# 1. UMAP of GSE171213 by cell type and condition
meta = pd.read_csv("data/processed/ml_input/cell_metadata.csv")
broad = {"T cells": "T cells", "Cytotoxic T cells": "T cells", "NK/cytotoxic lymphocytes": "NK / cytotoxic",
         "Proliferating cells": "Proliferating", "B cells": "B cells", "Plasma cells": "Plasma cells",
         "Neutrophils": "Neutrophils", "Monocytes/macrophages": "Mono/Macro", "Uncertain myeloid-like": "Mono/Macro",
         "Mast cells": "Mast cells", "Endothelial cells": "Endothelial", "Fibroblasts": "Fibroblasts",
         "Pericytes/smooth-muscle-like": "Pericytes/SMC", "Epithelial cells": "Epithelial", "Erythroid cells": "Erythroid"}
meta["label"] = meta.cell_type.map(broad)
pal = dict(zip(sorted(meta.label.unique()), plt.cm.tab20.colors))
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
for lab, d in meta.groupby("label"):
    ax[0].scatter(d.UMAP_1, d.UMAP_2, s=0.6, color=pal[lab], rasterized=True)
    ax[0].text(d.UMAP_1.median(), d.UMAP_2.median(), lab, fontsize=9, ha="center", weight="bold",
               color=INK, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.75))
idx = np.random.default_rng(0).permutation(len(meta))
m2 = meta.iloc[idx]
ax[1].scatter(m2.UMAP_1, m2.UMAP_2, s=0.6, c=np.where(m2.condition == "Healthy", HC, PD), rasterized=True)
for c, col in (("Healthy", HC), ("Periodontitis", PD)):
    ax[1].scatter([], [], color=col, s=30, label=f"{c} ({(meta.condition == c).sum():,} cells)")
ax[1].legend(frameon=False, loc="upper right")
for a, t in zip(ax, ("Annotated cell types", "Condition")):
    a.set_title(t, loc="left", weight="bold"); a.set_xticks([]); a.set_yticks([]); a.grid(False)
    a.set_xlabel("UMAP 1"); a.set_ylabel("UMAP 2")
save(fig, "fig1_umap")

# 2. composition: per-donor proportions for key types in both cohorts
cb = pd.read_csv(f"{T}/composition_by_donor_both_cohorts.csv")
ct_tests = pd.read_csv(f"{T}/composition_two_scRNA_cohorts.csv")
types = ["Plasma cells", "Monocytes/macrophages", "Endothelial cells", "T/NK cells",
         "B cells", "Pericytes/smooth muscle", "Epithelial cells", "Fibroblasts"]
fig, axes = plt.subplots(2, 8, figsize=(17, 6.2), sharey=False)
rng = np.random.default_rng(1)
for r, ds in enumerate(("GSE171213", "GSE164241")):
    for c, ct in enumerate(types):
        a = axes[r, c]; d = cb[(cb.dataset == ds) & (cb.broad == ct)]
        for i, (cond, col) in enumerate((("Healthy", HC), ("Periodontitis", PD))):
            v = 100 * d[d.condition == cond].prop
            a.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, color=col, s=26, edgecolor="white", lw=0.6, zorder=3)
            a.hlines(v.median(), i - 0.25, i + 0.25, color=INK, lw=2, zorder=4)
        row = ct_tests[(ct_tests.dataset == ds) & (ct_tests.cell_type == ct)].iloc[0]
        a.set_title(f"{ct.replace('Monocytes/macrophages','Mono/Macro').replace('Pericytes/smooth muscle','Pericytes/SMC')}\nFDR={row.FDR_BH:.2f}",
                    fontsize=9.5, weight="bold" if row.FDR_BH < 0.05 else "normal")
        a.set_xticks([0, 1]); a.set_xticklabels(["HC", "PD"]); a.set_xlim(-0.6, 1.6)
        if c == 0: a.set_ylabel(f"{ds}\n% of cells", weight="bold")
fig.suptitle("Cell-type proportions per donor (Mann-Whitney, BH-adjusted); bold = FDR < 0.05", x=0.01, ha="left", weight="bold")
fig.tight_layout(); save(fig, "fig2_composition_two_cohorts")

# 3. DE counts: paper-style Wilcoxon (cells as replicates) vs donor pseudobulk
w = pd.read_csv(f"{T}/DE_paper_style_wilcoxon.csv"); pb = pd.read_csv(f"{T}/DE_all_cell_types.csv")
wc = w[w.p_val_adj < 0.05].groupby("cell_type").size(); pc = pb[pb.FDR < 0.05].groupby("cell_type").size()
dd = pd.DataFrame({"Wilcoxon (cell-level, paper method)": wc, "edgeR pseudobulk (donor-level)": pc}).fillna(0)
dd = dd.sort_values("Wilcoxon (cell-level, paper method)")
fig, ax = plt.subplots(figsize=(9, 6))
y = np.arange(len(dd))
ax.barh(y + 0.2, dd.iloc[:, 0], height=0.38, color="#9b9a95", label=dd.columns[0])
ax.barh(y - 0.2, dd.iloc[:, 1], height=0.38, color=PD, label=dd.columns[1])
for i, (a_, b_) in enumerate(dd.values):
    ax.text(a_ + 8, i + 0.2, f"{int(a_)}", va="center", fontsize=9, color=INK2)
    ax.text(b_ + 8, i - 0.2, f"{int(b_)}", va="center", fontsize=9, color=INK)
ax.set_yticks(y); ax.set_yticklabels(dd.index); ax.set_xlabel("Genes significant (adjusted p < 0.05)")
ax.legend(frameon=False, loc="lower right"); ax.grid(axis="y", visible=False)
ax.set_title(f"Treating cells as replicates inflates DE: {int(dd.iloc[:,0].sum()):,} vs {int(dd.iloc[:,1].sum())} gene hits", loc="left", weight="bold")
save(fig, "fig3_de_method_comparison")

# 4. candidate gene replication heatmap
rep = pd.read_csv(f"{T}/candidate_genes_replication.csv")
lab = rep.gene + "  (" + rep.cell_type.str.replace(" cells", "") + ")"
M = np.column_stack([np.sign(rep.logFC_GSE171213) * np.minimum(abs(rep.logFC_GSE171213), 3),
                     np.sign(rep.logFC_GSE164241) * np.minimum(abs(rep.logFC_GSE164241), 3),
                     np.sign(rep.log2FC_bulk_GSE10334) * np.minimum(abs(rep.log2FC_bulk_GSE10334), 3)])
P = np.column_stack([rep.FDR_GSE171213, rep.p_GSE164241, rep.FDR_bulk])
cmap = LinearSegmentedColormap.from_list("div", [HC, "#f0efec", PD])
fig, ax = plt.subplots(figsize=(7.4, 7.6))
im = ax.imshow(M, cmap=cmap, vmin=-3, vmax=3, aspect="auto")
for i in range(M.shape[0]):
    for j in range(3):
        if np.isnan(M[i, j]): ax.text(j, i, "n/a", ha="center", va="center", color=INK2, fontsize=8); continue
        s = "**" if P[i, j] < 0.01 else "*" if P[i, j] < 0.05 else ""
        ax.text(j, i, s, ha="center", va="center", color=INK, fontsize=13, weight="bold")
ax.set_xticks(range(3)); ax.set_xticklabels(["GSE171213\n(discovery)", "GSE164241\n(scRNA)", "GSE10334\n(bulk, 247)"])
ax.set_yticks(range(len(lab))); ax.set_yticklabels(lab); ax.grid(False)
robust = rep.agree_GSE164241 & rep.agree_bulk & (rep.p_GSE164241 < 0.1) & (rep.FDR_bulk < 0.05)
cb_ = fig.colorbar(im, ax=ax, shrink=0.6)
ax.set_yticklabels([l + ("  [replicated]" if r else "") for l, r in zip(lab, robust)])
for t, r in zip(ax.get_yticklabels(), robust):
    t.set_fontweight("bold" if r else "normal"); t.set_color(INK if r else INK2); cb_.set_label("log2 fold change PD vs HC (capped at ±3)")
ax.set_title("Replication of the 14 pseudobulk hits\n* p/FDR<0.05  ** <0.01; bold = replicated in both", loc="left", weight="bold", fontsize=11)
save(fig, "fig4_candidate_replication")

# 5. enrichment dot plot (selected, ambient-robust terms)
en = pd.read_csv(f"{T}/enrichment_GO_KEGG_by_cell_type.csv")
pick = [("Neutrophils", "up", "Osteoclast differentiation"), ("Neutrophils", "up", "NF-kappa B signaling pathway"),
        ("Neutrophils", "up", "Response To Cytokine (GO:0034097)"), ("Neutrophils", "down", "Defense Response To Symbiont (GO:0140546)"),
        ("Plasma cells", "up", "Protein processing in endoplasmic reticulum"), ("Plasma cells", "down", "Collagen Fibril Organization (GO:0030199)"),
        ("Monocytes/macrophages", "up", "Phagosome"), ("Monocytes/macrophages", "up", "Response To Endoplasmic Reticulum Stress (GO:0034976)"),
        ("Mast cells", "up", "Regulation Of I-kappaB kinase/NF-kappaB Signaling (GO:0043122)"),
        ("Mast cells", "up", "Regulation Of Tumor Necrosis Factor Production (GO:0032680)"),
        ("Fibroblasts", "up", "Aerobic Electron Transport Chain (GO:0019646)"), ("Fibroblasts", "down", "Extracellular Matrix Organization (GO:0030198)"),
        ("Endothelial cells", "down", "Prostaglandin Metabolic Process (GO:0006693)"), ("Endothelial cells", "up", "Regulation Of Apoptotic Process (GO:0042981)"),
        ("NK/cytotoxic lymphocytes", "up", "Cellular Response To Tumor Necrosis Factor (GO:0071356)"),
        ("T cells", "down", "Response To Tumor Necrosis Factor (GO:0034612)")]
rows = []
for ct, dr, term in pick:
    r = en[(en.cell_type == ct) & (en.direction == dr) & (en.Term == term)]
    if len(r): rows.append(r.iloc[0])
e = pd.DataFrame(rows)
e["label"] = e.Term.str.replace(r" \(GO:\d+\)", "", regex=True) + "  —  " + e.cell_type.str.replace("Monocytes/macrophages", "Mono/Macro").str.replace("NK/cytotoxic lymphocytes", "NK") + " " + np.where(e.direction == "up", "(up)", "(down)")
e["k"] = e.Overlap.str.split("/").str[0].astype(int); e["mlp"] = -np.log10(e["Adjusted P-value"])
e = e.iloc[::-1]
fig, ax = plt.subplots(figsize=(10, 7))
ax.scatter(e.mlp, range(len(e)), s=e.k * 14, c=np.where(e.direction == "up", PD, HC), edgecolor="white", lw=1, zorder=3)
ax.set_yticks(range(len(e))); ax.set_yticklabels(e.label, fontsize=9.5); ax.grid(axis="y", visible=False)
ax.set_xlabel("−log10 adjusted p (Enrichr)")
for k in (5, 10, 20): ax.scatter([], [], s=k * 14, color="#bbb", label=f"{k} genes")
ax.scatter([], [], color=PD, s=60, label="up in PD"); ax.scatter([], [], color=HC, s=60, label="down in PD")
ax.legend(frameon=False, loc="lower right", fontsize=9, labelspacing=1.1)
ax.set_title("Selected GO / KEGG terms enriched among cell-type DE genes", loc="left", weight="bold")
save(fig, "fig5_enrichment")

# 6. bulk validation: volcano + ML comparison
bulk = pd.read_csv(f"{T}/DE_bulk_GSE10334.csv")
fig, ax = plt.subplots(figsize=(7.5, 6))
sig = (bulk.FDR < 0.05) & (bulk.log2FC.abs() > 1)
ax.scatter(bulk.log2FC[~sig], -np.log10(bulk.p[~sig]), s=3, color="#c9c8c3", rasterized=True)
ax.scatter(bulk.log2FC[sig], -np.log10(bulk.p[sig]), s=8, c=np.where(bulk.log2FC[sig] > 0, PD, HC))
for _, r in bulk[sig].nsmallest(14, "p").iterrows():
    ax.annotate(r.gene, (r.log2FC, -np.log10(r.p)), fontsize=8.5, xytext=(4, 2), textcoords="offset points")
for g in ("SELL", "CYR61", "CXCL6", "THBS1"):
    r = bulk[bulk.gene == g].iloc[0]
    ax.annotate(g, (r.log2FC, -np.log10(r.p)), fontsize=9, weight="bold", color=INK,
                xytext=(8, -12), textcoords="offset points", arrowprops=dict(arrowstyle="-", color=INK2))
ax.set_xlabel("log2 fold change (diseased vs healthy site)"); ax.set_ylabel("−log10 p (Welch t-test)")
ax.set_title(f"GSE10334 bulk biopsies: {int((bulk.FDR<0.05).sum()):,} genes FDR<0.05, {int(sig.sum())} with |log2FC|>1",
             loc="left", weight="bold", fontsize=11)
save(fig, "fig6_bulk_volcano")

cv = pd.read_csv(f"{T}/ml_bulk_cv_summary.csv", index_col=0); ext = pd.read_csv(f"{T}/ml_bulk_external_test.csv", index_col=0)
names = {"Majority class (always periodontitis)": "Majority class",
         "Signature score, no training (scRNA-derived)": "scRNA signature score",
         "Best single gene (t-test in training folds)": "Best single gene",
         "Bulk DE score: mean z(top-25 up) - z(top-25 down)": "DE gene score (50 genes)",
         "Logistic regression, scRNA signature genes": "LR, scRNA genes",
         "Logistic regression, top-50 genes": "LR, top-50 genes",
         "Linear SVM, top-50 genes": "Linear SVM, top-50",
         "Random forest, top-50 genes": "Random forest, top-50"}
group = ["base", "hyp", "hyp", "hyp", "ml", "ml", "ml", "ml"]
gcol = {"base": "#c9c8c3", "hyp": "#1baf7a", "ml": "#6250d6"}
fig, ax = plt.subplots(1, 2, figsize=(14, 5.6), sharey=True)
yy = np.arange(len(cv))[::-1]
for a, metric, title in ((ax[0], "ROC_AUC", "ROC-AUC"), (ax[1], "balanced_accuracy", "Balanced accuracy")):
    a.barh(yy, cv[f"{metric}_mean"], xerr=cv[f"{metric}_std"], color=[gcol[g] for g in group], height=0.6,
           error_kw=dict(ecolor=INK2, lw=1, capsize=3))
    exv = [ext.loc[m, metric] if m in ext.index else np.nan for m in cv.index]
    a.scatter(exv, yy, marker="D", color=INK, s=36, zorder=4, label="External test (GSE16134, 29 new patients)")
    for yv, v in zip(yy, cv[f"{metric}_mean"]): a.text(0.02, yv, f"{v:.2f}", va="center", color="white", weight="bold", fontsize=9.5)
    a.axvline(0.5, color=INK2, ls="--", lw=1); a.set_xlim(0, 1.02); a.set_title(title + " (patient-grouped 5-fold CV ×5)", loc="left", weight="bold")
    a.grid(axis="y", visible=False)
ax[0].set_yticks(yy); ax[0].set_yticklabels([names[m] for m in cv.index])
from matplotlib.patches import Patch
ax[1].legend(handles=[Patch(color=gcol["hyp"], label="Hypothesis-testing based (no/1-parameter model)"),
                      Patch(color=gcol["ml"], label="Machine learning"), Patch(color=gcol["base"], label="Baseline"),
                      plt.Line2D([], [], marker="D", color=INK, ls="", label="External test (29 new patients)")],
             frameon=False, loc="upper center", bbox_to_anchor=(-0.05, -0.1), ncol=4, fontsize=10)
save(fig, "fig7_ml_bulk_comparison")

# ROC curves (CV out-of-fold, repeat 0)
pk = pd.read_pickle("data/external/ml_plot_data.pkl")
fig, ax = plt.subplots(figsize=(6, 5.6))
for (m, (fpr, tpr, _)), col in zip([(k, v) for k, v in pk["roc_cv"].items() if "Majority" not in k],
                                   ["#9b9a95", "#1baf7a", "#008300", CAT[0], CAT[6], CAT[1], CAT[4]]):
    ax.plot(fpr, tpr, lw=2, color=col, label=names[m])
ax.plot([0, 1], [0, 1], ls="--", color=INK2, lw=1)
ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate"); ax.legend(frameon=False, fontsize=8.5, loc="lower right")
ax.set_title("GSE10334 out-of-fold ROC curves", loc="left", weight="bold")
save(fig, "fig8_roc_bulk")

# 7. leakage: single-cell ML
s = json.load(open(f"{T}/cross_dataset_summary.json"))["scRNA_ml"]
fig, ax = plt.subplots(figsize=(7.5, 4.4))
labels = ["Random cell split\n(same donors in train & test)", "Leave-one-donor-out\n(cell level)", "Leave-one-donor-out\n(donor vote, n=9)"]
vals = [s["random_cell_split_AUC"], s["LODO_cell_AUC"], s["LODO_donor_AUC"]]
ax.bar(range(3), vals, color=["#e34948", "#6250d6", "#6250d6"], width=0.55)
for i, v in enumerate(vals): ax.text(i, v + 0.02, f"AUC {v:.2f}", ha="center", weight="bold")
ax.axhline(0.5, ls="--", color=INK2, lw=1); ax.set_ylim(0, 1.08); ax.set_xticks(range(3)); ax.set_xticklabels(labels, fontsize=9.5)
ax.set_ylabel("ROC-AUC"); ax.grid(axis="x", visible=False)
ax.set_title("Cell-level logistic regression, GSE171213: donor leakage", loc="left", weight="bold")
save(fig, "fig9_leakage")
print("figures written")
