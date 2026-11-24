"""(a) Compare GSE171213 with the independent scRNA cohort GSE164241 and the bulk cohort
GSE10334: composition, candidate genes and genome-wide tissue-level fold changes.
(b) Single-cell ML: show why cells must be split by donor (leakage), and what a
donor-level classifier achieves with only 9 donors."""
import json
import numpy as np
import pandas as pd
from scipy import stats, io
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, LeaveOneGroupOut
from sklearn.metrics import roc_auc_score, accuracy_score, balanced_accuracy_score

tab = "results/tables"
broad = {"T cells": "T/NK cells", "Cytotoxic T cells": "T/NK cells", "NK/cytotoxic lymphocytes": "T/NK cells",
         "Proliferating cells": "T/NK cells", "B cells": "B cells", "Plasma cells": "Plasma cells",
         "Neutrophils": "Neutrophils", "Monocytes/macrophages": "Monocytes/macrophages",
         "Uncertain myeloid-like": "Monocytes/macrophages", "Mast cells": "Mast cells",
         "Endothelial cells": "Endothelial cells", "Fibroblasts": "Fibroblasts",
         "Pericytes/smooth-muscle-like": "Pericytes/smooth muscle", "Epithelial cells": "Epithelial cells",
         "Erythroid cells": None}

# ---------------- composition in both scRNA cohorts at the broad level
meta = pd.read_csv("data/processed/ml_input/cell_metadata.csv")
meta["broad"] = meta.cell_type.map(broad)
c1 = meta.dropna(subset=["broad"]).groupby(["donor", "condition", "broad"]).size().unstack(fill_value=0)
c1 = c1.div(c1.sum(1), axis=0).stack().rename("prop").reset_index()
c2 = pd.read_csv(f"{tab}/GSE164241_composition_by_donor.csv").rename(columns={"cell_type": "broad"})
rows = []
for name, c in (("GSE171213", c1), ("GSE164241", c2)):
    for ct, d in c.groupby("broad"):
        h, p = d[d.condition == "Healthy"].prop, d[d.condition == "Periodontitis"].prop
        rows.append({"dataset": name, "cell_type": ct, "healthy_donors": len(h), "pd_donors": len(p),
                     "healthy_median_pct": 100 * h.median(), "pd_median_pct": 100 * p.median(),
                     "log2_ratio": np.log2((p.median() + 1e-3) / (h.median() + 1e-3)),
                     "p_value": stats.mannwhitneyu(h, p).pvalue})
comp = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
comp["FDR_BH"] = comp.groupby("dataset").p_value.transform(lambda p: multipletests(p, method="fdr_bh")[1])
comp.to_csv(f"{tab}/composition_two_scRNA_cohorts.csv", index=False)
pd.concat([c1.assign(dataset="GSE171213"), c2[["donor", "condition", "broad", "prop"]].assign(dataset="GSE164241")]) \
    .to_csv(f"{tab}/composition_by_donor_both_cohorts.csv", index=False)
wide = comp.pivot(index="cell_type", columns="dataset", values=["log2_ratio", "FDR_BH"])
print(wide.round(3))

# ---------------- candidate genes (GSE171213 cell-type pseudobulk FDR<0.05) in GSE164241
cand = pd.read_csv(f"{tab}/DE_all_cell_types.csv").query("FDR < 0.05")
cand["broad"] = cand.cell_type.map(broad)
v2 = pd.read_csv(f"{tab}/GSE164241_DE_celltype_pseudobulk.csv")
bulkcand = pd.read_csv(f"{tab}/candidate_genes_bulk_validation.csv")
out = cand.merge(v2[["gene", "cell_type", "logFC", "PValue", "FDR"]].rename(
    columns={"cell_type": "broad", "logFC": "logFC_GSE164241", "PValue": "p_GSE164241", "FDR": "FDR_GSE164241"}),
    on=["gene", "broad"], how="left")
out = out.merge(bulkcand[["gene", "cell_type", "log2FC", "FDR_bulk"]], on=["gene", "cell_type"], how="left")
out["agree_GSE164241"] = np.sign(out.logFC) == np.sign(out.logFC_GSE164241)
out["agree_bulk"] = np.sign(out.logFC) == np.sign(out.log2FC)
out = out[["cell_type", "gene", "logFC", "FDR", "logFC_GSE164241", "p_GSE164241", "agree_GSE164241",
           "log2FC", "FDR_bulk", "agree_bulk"]].rename(columns={"logFC": "logFC_GSE171213", "FDR": "FDR_GSE171213",
                                                               "log2FC": "log2FC_bulk_GSE10334"})
out.to_csv(f"{tab}/candidate_genes_replication.csv", index=False)
print(out.round(3).to_string())

# ---------------- genome-wide tissue-level fold-change concordance
t1 = pd.read_csv(f"{tab}/DE_tissue_pseudobulk.csv")[["gene", "logFC", "FDR"]]
t2 = pd.read_csv(f"{tab}/GSE164241_DE_tissue_pseudobulk.csv")[["gene", "logFC", "FDR"]]
b = pd.read_csv(f"{tab}/DE_bulk_GSE10334.csv")[["gene", "log2FC", "FDR"]]
m = t1.merge(t2, on="gene", suffixes=("_A", "_B")).merge(b.rename(columns={"log2FC": "logFC_bulk", "FDR": "FDR_bulk"}), on="gene")
m.to_csv(f"{tab}/tissue_logFC_three_cohorts.csv", index=False)
conc = {}
for a, bb in (("A", "B"), ("A", "bulk"), ("B", "bulk")):
    sel = m[(m[f"FDR_{a}"] < 0.05) & (m[f"FDR_{bb}"] < 0.05)]
    conc[f"{a}_vs_{bb}"] = {"rho_all": stats.spearmanr(m[f"logFC_{a}"], m[f"logFC_{bb}"]).statistic,
                            "n_both_sig": len(sel),
                            "direction_agreement_both_sig": float((np.sign(sel[f"logFC_{a}"]) == np.sign(sel[f"logFC_{bb}"])).mean())}
# replicated genes: significant & same direction in all three
rep = m[(m.FDR_A < 0.05) & (m.FDR_B < 0.05) & (m.FDR_bulk < 0.05) &
        (np.sign(m.logFC_A) == np.sign(m.logFC_B)) & (np.sign(m.logFC_A) == np.sign(m.logFC_bulk))]
rep = rep.assign(mean_logFC=rep[["logFC_A", "logFC_B", "logFC_bulk"]].mean(1)).sort_values("mean_logFC")
rep.to_csv(f"{tab}/replicated_genes_three_cohorts.csv", index=False)
conc["replicated_up"] = int((rep.logFC_A > 0).sum()); conc["replicated_down"] = int((rep.logFC_A < 0).sum())
print(json.dumps(conc, indent=1))
print("Replicated UP:", rep[rep.logFC_A > 0].sort_values("logFC_bulk", ascending=False).gene.head(30).tolist())
print("Replicated DOWN:", rep[rep.logFC_A < 0].sort_values("logFC_bulk").gene.head(20).tolist())

# ---------------- single-cell ML
X = io.mmread("data/processed/ml_input/lognorm_hvg_cells_x_genes.mtx").tocsr()
y = (meta.condition == "Periodontitis").astype(int).values
donor = meta.donor.values
def lr():
    return LogisticRegression(C=0.05, max_iter=3000, class_weight="balanced")
res = {}
# (1) random cell split: cells of the same donor in train and test -> leakage
aucs = []
for tr, te in StratifiedKFold(5, shuffle=True, random_state=0).split(X, y):
    p = lr().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]; aucs.append(roc_auc_score(y[te], p))
res["random_cell_split_AUC"] = float(np.mean(aucs))
# (2) leave-one-donor-out: the held-out donor is never seen in training
donor_rows, cell_p = [], np.zeros(len(y))
for tr, te in LeaveOneGroupOut().split(X, y, donor):
    p = lr().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]; cell_p[te] = p
    donor_rows.append({"donor": donor[te][0], "true": int(y[te][0]), "mean_prob_PD": float(p.mean()),
                       "frac_cells_called_PD": float((p > 0.5).mean())})
dd = pd.DataFrame(donor_rows)
res["LODO_cell_AUC"] = float(roc_auc_score(y, cell_p))
res["LODO_cell_accuracy"] = float(accuracy_score(y, cell_p > 0.5))
res["LODO_donor_accuracy"] = float(accuracy_score(dd.true, dd.mean_prob_PD > 0.5))
res["LODO_donor_AUC"] = float(roc_auc_score(dd.true, dd.mean_prob_PD))
dd.to_csv(f"{tab}/ml_scRNA_LODO_donor_predictions.csv", index=False)
# (3) donor-level pseudobulk, 9 samples, leave-one-out
pb = pd.read_csv("data/processed/ml_input/pseudobulk_donor_logcpm.csv", index_col=0)
Xd = pb.T.values; yd = np.array([1 if d.startswith("PD") else 0 for d in pb.columns])
pp = np.zeros(len(yd))
for i in range(len(yd)):
    tr = np.arange(len(yd)) != i
    t_, _ = stats.ttest_ind(Xd[tr][yd[tr] == 1], Xd[tr][yd[tr] == 0], equal_var=False)
    idx = np.argsort(-np.abs(np.nan_to_num(t_)))[:50]
    mu, sd = Xd[tr][:, idx].mean(0), Xd[tr][:, idx].std(0) + 1e-9
    m_ = LogisticRegression(C=0.1, max_iter=3000).fit((Xd[tr][:, idx] - mu) / sd, yd[tr])
    pp[i] = m_.predict_proba(((Xd[i, idx] - mu) / sd)[None])[0, 1]
res["pseudobulk_LOO_accuracy"] = float(accuracy_score(yd, pp > 0.5))
res["pseudobulk_LOO_AUC"] = float(roc_auc_score(yd, pp))
print(json.dumps(res, indent=1)); print(dd.round(3).to_string())
json.dump({"concordance": conc, "scRNA_ml": res}, open(f"{tab}/cross_dataset_summary.json", "w"), indent=1)
