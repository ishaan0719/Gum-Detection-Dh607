"""Parse the bulk gingival-tissue microarray cohorts (GSE10334, GSE16134; Affymetrix
HG-U133 Plus 2.0), collapse probes to gene symbols and build an external cohort that
does not share samples with GSE10334."""
import gzip, re
from pathlib import Path
import numpy as np
import pandas as pd

ext = Path("data/external")

def read_series(gse):
    meta, rows = {}, []
    with gzip.open(ext / f"{gse}_series_matrix.txt.gz", "rt") as fh:
        in_table = False
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True; continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                rows.append(line.rstrip("\n").split("\t")); continue
            if line.startswith("!Sample_"):
                parts = line.rstrip("\n").split("\t")
                key = parts[0]
                vals = [p.strip('"') for p in parts[1:]]
                meta.setdefault(key, []).append(vals)
    header = [h.strip('"') for h in rows[0]]
    expr = pd.DataFrame([r[1:] for r in rows[1:]], index=[r[0].strip('"') for r in rows[1:]],
                        columns=header[1:]).apply(pd.to_numeric, errors="coerce")
    titles = meta["!Sample_title"][0]
    pheno = [v for v in meta["!Sample_characteristics_ch1"] if "phenotype" in v[0].lower()][0]
    s = pd.DataFrame({"gsm": meta["!Sample_geo_accession"][0], "title": titles,
                      "phenotype": pheno})
    s["patient"] = s.title.str.extract(r"patient (\d+)", flags=re.I)[0]
    s["condition"] = np.where(s.phenotype.str.contains("Diseased", case=False),
                              "Periodontitis", "Healthy")
    s["cohort"] = gse
    return expr, s.set_index("gsm")

# The .annot file has a preamble; locate the real header.
if True:
    with gzip.open(ext / "GPL570.annot.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if line.startswith("ID\t"):
                skip = i; break
    annot = pd.read_csv(ext / "GPL570.annot.gz", sep="\t", skiprows=skip, low_memory=False)
probe2gene = annot.set_index("ID")["Gene symbol"].dropna()
probe2gene = probe2gene[~probe2gene.str.contains("///")]

def collapse(expr):
    e = expr.loc[expr.index.intersection(probe2gene.index)]
    if e.max().max() > 100:  # MAS5-style scale -> log2
        e = np.log2(e.clip(lower=1))
    e = e.assign(gene=probe2gene.loc[e.index].values)
    # keep the most variable probe per gene
    e["var"] = e.drop(columns="gene").var(axis=1)
    e = e.sort_values("var", ascending=False).drop_duplicates("gene")
    return e.drop(columns="var").set_index("gene")

out = {}
for gse in ("GSE10334", "GSE16134"):
    expr, s = read_series(gse)
    g = collapse(expr)
    out[gse] = (g, s)
    print(gse, "samples", s.shape[0], "genes", g.shape[0],
          s.condition.value_counts().to_dict(), "patients", s.patient.nunique(),
          "max expr", round(float(g.max().max()), 1))

g10, s10 = out["GSE10334"]; g16, s16 = out["GSE16134"]
# GSE16134 re-deposits most GSE10334 arrays under new GSM IDs (different processing),
# so GSM matching finds nothing. Detect duplicates from within-cohort standardised,
# rank-transformed profiles and drop every GSE16134 patient with a duplicated array.
genes = g10.index.intersection(g16.index)
r10, r16 = g10.loc[genes].rank(), g16.loc[genes].rank()
zs = lambda x: x.sub(x.mean(1), axis=0).div(x.std(1) + 1e-9, axis=0)
hv = r10.var(1) > r10.var(1).median()
A, B = zs(r10)[hv], zs(r16)[hv]
A, B = (A - A.mean()) / A.std(), (B - B.mean()) / B.std()
best = ((B.T.values @ A.values) / A.shape[0]).max(1)
dup = pd.Series(best > 0.8, index=g16.columns)
pd.DataFrame({"best_corr_to_GSE10334": best, "duplicate": dup.values},
             index=g16.columns).to_csv(ext / "GSE16134_duplicate_check.csv")
dup_patients = set(s16.loc[dup[dup].index, "patient"])
ext_keep = s16.index[~s16.patient.isin(dup_patients)]
print("GSE16134 arrays duplicating GSE10334:", int(dup.sum()))
print("Independent GSE16134 samples:", len(ext_keep),
      s16.loc[ext_keep, "condition"].value_counts().to_dict(),
      "patients", s16.loc[ext_keep, "patient"].nunique())

g10.to_csv(ext / "GSE10334_gene_expr.csv.gz")
s10.to_csv(ext / "GSE10334_samples.csv")
g16[ext_keep].to_csv(ext / "GSE16134_independent_gene_expr.csv.gz")
s16.loc[ext_keep].to_csv(ext / "GSE16134_independent_samples.csv")
