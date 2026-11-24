"""Independent validation in bulk gingival biopsies and the ML comparison.

Step 1 (hypothesis testing first): do genes / signatures found in the single-cell data
differ between diseased and healthy sites in GSE10334 (247 biopsies, 90 patients)?
Step 2 (ML): can LR / linear SVM (and RF for reference) classify sites, using
patient-grouped cross-validation, and do they beat the training-free signature score?
Step 3: external test on GSE16134 patients not present in GSE10334.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, precision_score,
                             recall_score, f1_score, roc_auc_score, roc_curve)
import gseapy as gp

rng = np.random.default_rng(0)
ext = Path("data/external"); tab = Path("results/tables"); tab.mkdir(exist_ok=True, parents=True)

X10 = pd.read_csv(ext / "GSE10334_gene_expr.csv.gz", index_col=0)
s10 = pd.read_csv(ext / "GSE10334_samples.csv", index_col=0).loc[X10.columns]
X16 = pd.read_csv(ext / "GSE16134_independent_gene_expr.csv.gz", index_col=0)
s16 = pd.read_csv(ext / "GSE16134_independent_samples.csv", index_col=0).loc[X16.columns]
y10 = (s10.condition == "Periodontitis").astype(int).values
y16 = (s16.condition == "Periodontitis").astype(int).values
g10 = s10.patient.values

# keep reasonably expressed genes
expressed = X10.index[X10.mean(1) > X10.mean(1).quantile(0.25)]
X10, X16 = X10.loc[expressed], X16.loc[X16.index.intersection(expressed)]

# ---------------------------------------------------------------- 1. bulk DE
d, h = X10.loc[:, y10 == 1], X10.loc[:, y10 == 0]
t, p = stats.ttest_ind(d, h, axis=1, equal_var=False)
bulk = pd.DataFrame({"gene": X10.index, "log2FC": d.mean(1) - h.mean(1), "t": t, "p": p})
bulk["FDR"] = multipletests(bulk.p, method="fdr_bh")[1]
bulk = bulk.sort_values("p").reset_index(drop=True)
bulk.to_csv(tab / "DE_bulk_GSE10334.csv", index=False)
n_sig = int(((bulk.FDR < 0.05) & (bulk.log2FC.abs() > 1)).sum())
print("Bulk DE: FDR<0.05:", int((bulk.FDR < 0.05).sum()), " FDR<0.05 & |log2FC|>1:", n_sig)

# patient-paired check: patients with both a healthy and a diseased site
paired = s10.groupby("patient").condition.nunique()
paired_pat = paired.index[paired == 2]
print("Patients with both healthy and diseased sites:", len(paired_pat))

# --------------------------------------------- 2. single-cell signatures in bulk
tissue = pd.read_csv(tab / "DE_tissue_pseudobulk.csv")
tissue = tissue[tissue.gene.isin(X10.index)]
sig_t = tissue[tissue.FDR < 0.05]
signatures = {
    "scRNA tissue pseudobulk (top 50 up / 50 down)": (
        sig_t[sig_t.logFC > 0].nsmallest(50, "PValue").gene.tolist(),
        sig_t[sig_t.logFC < 0].nsmallest(50, "PValue").gene.tolist()),
}
ct = pd.read_csv(tab / "DE_all_cell_types.csv")
ct = ct[(ct.FDR < 0.05) & ct.gene.isin(X10.index)]
signatures["scRNA cell-type pseudobulk hits (FDR<0.05)"] = (
    ct[ct.logFC > 0].gene.tolist(), ct[ct.logFC < 0].gene.tolist())
w = pd.read_csv(tab / "DE_paper_style_wilcoxon.csv")
w = w[(w.p_val_adj < 0.05) & w.gene.isin(X10.index) & ~w.gene.str.match(r"^(RPL|RPS|MT-)")]
cons = w.groupby("gene").avg_log2FC.agg(["count", lambda v: np.sign(v).mean()])
cons.columns = ["n", "dir"]
signatures["Paper-style Wilcoxon consensus (>=3 cell types)"] = (
    cons[(cons.n >= 3) & (cons.dir == 1)].index.tolist(),
    cons[(cons.n >= 3) & (cons.dir == -1)].index.tolist())

def zscore(X):
    return X.sub(X.mean(1), axis=0).div(X.std(1) + 1e-9, axis=0)

def score(X, up, down):
    Z = zscore(X)
    up = [g for g in up if g in Z.index]; down = [g for g in down if g in Z.index]
    s = 0
    if up: s = s + Z.loc[up].mean()
    if down: s = s - Z.loc[down].mean()
    return s

def boot_auc(y, s, n=2000):
    vals = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i])) == 2:
            vals.append(roc_auc_score(y[i], s[i]))
    return np.percentile(vals, [2.5, 97.5])

sig_rows, sig_scores = [], {}
for name, (up, down) in signatures.items():
    s_10 = score(X10, up, down).values; s_16 = score(X16, up, down).values
    sig_scores[name] = (s_10, s_16)
    pat = pd.DataFrame({"s": s_10, "pat": g10, "c": y10}).groupby(["pat", "c"]).s.mean().unstack()
    pat = pat.loc[paired_pat]
    lo, hi = boot_auc(y10, s_10)
    lo16, hi16 = boot_auc(y16, s_16)
    sig_rows.append({
        "signature": name, "n_up": len([g for g in up if g in X10.index]),
        "n_down": len([g for g in down if g in X10.index]),
        "GSE10334_AUC": roc_auc_score(y10, s_10), "GSE10334_AUC_CI": f"{lo:.2f}-{hi:.2f}",
        "GSE10334_MannWhitney_p": stats.mannwhitneyu(s_10[y10 == 1], s_10[y10 == 0]).pvalue,
        "paired_patients": len(pat),
        "paired_frac_diseased_higher": float((pat[1] > pat[0]).mean()),
        "paired_Wilcoxon_p": stats.wilcoxon(pat[1], pat[0]).pvalue,
        "GSE16134_indep_AUC": roc_auc_score(y16, s_16), "GSE16134_AUC_CI": f"{lo16:.2f}-{hi16:.2f}",
    })
sig_df = pd.DataFrame(sig_rows)
sig_df.to_csv(tab / "signature_validation_bulk.csv", index=False)
print(sig_df.to_string())

# per-gene validation of the 14 cell-type pseudobulk hits
cand = pd.read_csv(tab / "DE_all_cell_types.csv").query("FDR < 0.05")
cand = cand.merge(bulk, on="gene", how="left", suffixes=("_scRNA", "_bulk"))
cand["direction_agrees"] = np.sign(cand.logFC) == np.sign(cand.log2FC)
cand[["cell_type", "gene", "logFC", "FDR_scRNA", "log2FC", "p", "FDR_bulk", "direction_agrees"]] \
    .to_csv(tab / "candidate_genes_bulk_validation.csv", index=False)
print(cand[["cell_type", "gene", "logFC", "log2FC", "FDR_bulk", "direction_agrees"]].to_string())

# genome-wide concordance between scRNA tissue logFC and bulk log2FC
conc = tissue.merge(bulk, on="gene")
rho_all = stats.spearmanr(conc.logFC, conc.log2FC)
top = conc[conc.FDR_x < 0.05] if "FDR_x" in conc else conc[conc.FDR < 0.05]
rho_sig = stats.spearmanr(top.logFC, top.log2FC)
agree_sig = float((np.sign(top.logFC) == np.sign(top.log2FC)).mean())
conc[["gene", "logFC", "log2FC"]].rename(columns={"logFC": "scRNA_logFC", "log2FC": "bulk_log2FC"}) \
    .to_csv(tab / "scRNA_vs_bulk_logFC.csv", index=False)
print(f"Concordance all genes rho={rho_all.statistic:.2f}; scRNA FDR<0.05 genes rho={rho_sig.statistic:.2f}, "
      f"direction agreement={agree_sig:.2f} (n={len(top)})")

# enrichment of the bulk DE genes
up_b = bulk[(bulk.FDR < 0.05) & (bulk.log2FC > 1)].gene.tolist()
dn_b = bulk[(bulk.FDR < 0.05) & (bulk.log2FC < -1)].gene.tolist()
enr_rows = []
for direction, genes in (("up", up_b), ("down", dn_b)):
    if len(genes) >= 10:
        r = gp.enrichr(gene_list=genes, gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Human"],
                       organism="human", outdir=None).results
        r["direction"] = direction; enr_rows.append(r[r["Adjusted P-value"] < 0.05])
pd.concat(enr_rows).to_csv(tab / "enrichment_bulk_GSE10334.csv", index=False)

# ------------------------------------------------------------------- 3. ML
def metrics(y, prob, pred):
    return {"accuracy": accuracy_score(y, pred), "balanced_accuracy": balanced_accuracy_score(y, pred),
            "precision": precision_score(y, pred, zero_division=0), "recall": recall_score(y, pred),
            "F1": f1_score(y, pred), "ROC_AUC": roc_auc_score(y, prob)}

def select_top(Xtr, ytr, k):
    t, _ = stats.ttest_ind(Xtr[ytr == 1], Xtr[ytr == 0], axis=0, equal_var=False)
    return np.argsort(-np.abs(np.nan_to_num(t)))[:k]

genes = X10.index.intersection(X16.index)
A10 = X10.loc[genes].T.values; A16 = X16.loc[genes].T.values
sig_up, sig_dn = signatures["scRNA tissue pseudobulk (top 50 up / 50 down)"]
sig_idx = np.array([genes.get_loc(g) for g in sig_up + sig_dn if g in genes])

def std_fit(Xtr):
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    return lambda X: (X - mu) / sd

models = {
    "Majority class (always periodontitis)": None,
    "Signature score, no training (scRNA-derived)": "score",
    "Best single gene (t-test in training folds)": "best1",
    "Bulk DE score: mean z(top-25 up) - z(top-25 down)": "descore",
    "Logistic regression, scRNA signature genes": ("lr", "sig"),
    "Logistic regression, top-50 genes": ("lr", 50),
    "Linear SVM, top-50 genes": ("svm", 50),
    "Random forest, top-50 genes": ("rf", 50),
}

def make(kind):
    if kind == "lr":
        return LogisticRegression(C=0.1, class_weight="balanced", max_iter=5000)
    if kind == "svm":
        return LinearSVC(C=0.01, class_weight="balanced", max_iter=20000)
    return RandomForestClassifier(n_estimators=500, class_weight="balanced", random_state=0, n_jobs=-1)

def fit_predict(spec, Xtr, ytr, Xte, Xte_is_external=False):
    if spec is None:
        return np.full(len(Xte), 0.5) + rng.normal(0, 1e-6, len(Xte)), np.ones(len(Xte), int)
    if spec == "score":
        s_tr = score_arr(Xtr); s_te = score_arr(Xte)
        lr = LogisticRegression(class_weight="balanced").fit(s_tr[:, None], ytr)
        return s_te, lr.predict(s_te[:, None])
    if spec in ("best1", "descore"):
        t, _ = stats.ttest_ind(Xtr[ytr == 1], Xtr[ytr == 0], axis=0, equal_var=False)
        t = np.nan_to_num(t)
        if spec == "best1":
            j = np.argmax(np.abs(t)); f = lambda X: np.sign(t[j]) * X[:, j]
        else:
            u, dn = np.argsort(-t)[:25], np.argsort(t)[:25]
            mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
            f = lambda X: ((X - mu) / sd)[:, u].mean(1) - ((X - mu) / sd)[:, dn].mean(1)
        lr = LogisticRegression(class_weight="balanced").fit(f(Xtr)[:, None], ytr)
        return f(Xte), lr.predict(f(Xte)[:, None])
    kind, feat = spec
    idx = sig_idx if feat == "sig" else select_top(Xtr, ytr, feat)
    sc = std_fit(Xtr[:, idx])
    m = make(kind).fit(sc(Xtr[:, idx]), ytr)
    Xt = sc(Xte[:, idx])
    prob = m.decision_function(Xt) if kind == "svm" else m.predict_proba(Xt)[:, 1]
    return prob, m.predict(Xt)

up_i = np.array([genes.get_loc(g) for g in sig_up if g in genes])
dn_i = np.array([genes.get_loc(g) for g in sig_dn if g in genes])
def score_arr(X):  # z within the provided set of samples (no labels used)
    Z = (X - X.mean(0)) / (X.std(0) + 1e-9)
    return Z[:, up_i].mean(1) - Z[:, dn_i].mean(1)

cv_rows, roc_store = [], {}
for name, spec in models.items():
    for rep in range(5):
        cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=rep)
        oof_p = np.zeros(len(y10)); oof_y = np.zeros(len(y10), int)
        for tr, te in cv.split(A10, y10, g10):
            assert not set(g10[tr]) & set(g10[te])  # no patient in both train and test
            prob, pred = fit_predict(spec, A10[tr], y10[tr], A10[te])
            oof_p[te], oof_y[te] = prob, pred
            m = metrics(y10[te], prob, pred); m.update(model=name, repeat=rep)
            cv_rows.append(m)
        if rep == 0:
            roc_store[name] = roc_curve(y10, oof_p)
cv_df = pd.DataFrame(cv_rows)
summary = cv_df.groupby("model", sort=False).agg(["mean", "std"]).drop(columns="repeat")
summary.columns = [f"{a}_{b}" for a, b in summary.columns]
summary.to_csv(tab / "ml_bulk_cv_summary.csv")
print(summary.round(3)[[c for c in summary.columns if c.endswith("mean")]].to_string())

# Random sample-level CV (ignores patient) for the leakage comparison
from sklearn.model_selection import StratifiedKFold
leak = []
for rep in range(5):
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=rep).split(A10, y10):
        prob, pred = fit_predict(("lr", 50), A10[tr], y10[tr], A10[te])
        leak.append(roc_auc_score(y10[te], prob))
print("LR top-50 random-sample CV AUC:", np.mean(leak).round(3))

# external test: train on all of GSE10334, test on independent GSE16134 patients.
# Each cohort is gene-wise standardised separately (different array processing).
Z10 = (A10 - A10.mean(0)) / (A10.std(0) + 1e-9); Z16 = (A16 - A16.mean(0)) / (A16.std(0) + 1e-9)
ext_rows, ext_roc = [], {}
for name, spec in models.items():
    if spec is None:
        continue
    prob, pred = fit_predict(spec, Z10, y10, Z16)
    m = metrics(y16, prob, pred); lo, hi = boot_auc(y16, prob)
    m.update(model=name, AUC_95CI=f"{lo:.2f}-{hi:.2f}"); ext_rows.append(m)
    ext_roc[name] = roc_curve(y16, prob)
ext_df = pd.DataFrame(ext_rows).set_index("model")
ext_df.to_csv(tab / "ml_bulk_external_test.csv")
print(ext_df.round(3).to_string())

# LR coefficients on the full cohort (top-50 genes) for interpretation
idx = select_top(A10, y10, 50); sc = std_fit(A10[:, idx])
lr = make("lr").fit(sc(A10[:, idx]), y10)
coef = pd.DataFrame({"gene": genes[idx], "coef": lr.coef_[0]}).sort_values("coef")
coef.to_csv(tab / "ml_bulk_lr_top50_coefficients.csv", index=False)

pd.to_pickle({"roc_cv": roc_store, "roc_ext": ext_roc, "sig_scores": sig_scores,
              "y10": y10, "y16": y16, "patients10": g10}, ext / "ml_plot_data.pkl")
json.dump({"bulk_fdr05": int((bulk.FDR < 0.05).sum()), "bulk_fdr05_fc1": n_sig,
           "paired_patients": int(len(paired_pat)), "rho_all": rho_all.statistic,
           "rho_sig": rho_sig.statistic, "agree_sig": agree_sig, "n_sig_concord": int(len(top)),
           "leaky_auc": float(np.mean(leak)), "n_ext": int(len(y16)),
           "n_ext_healthy": int((y16 == 0).sum()), "n_ext_patients": int(s16.patient.nunique())},
          open(tab / "bulk_summary.json", "w"), indent=1)
