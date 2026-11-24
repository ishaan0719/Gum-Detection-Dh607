"""GO Biological Process / KEGG over-representation (Enrichr) for genes that are
differentially expressed between periodontitis and healthy tissue within each major
cell type (paper-style Wilcoxon, Bonferroni p_adj < 0.05, |log2FC| > 0.5)."""
from pathlib import Path
import time
import pandas as pd
import gseapy as gp

w = pd.read_csv("results/tables/DE_paper_style_wilcoxon.csv")
sig = w[(w.p_val_adj < 0.05) & (w.avg_log2FC.abs() > 0.5)]
# drop ribosomal/mito/immunoglobulin/haemoglobin genes that mostly reflect ambient RNA
noise = r"^(RPL|RPS|MT-|IGH|IGK|IGL|HBA|HBB|MIR|AC\d|AL\d|LINC)"
sig = sig[~sig.gene.str.match(noise)]

libs = ["GO_Biological_Process_2023", "KEGG_2021_Human"]
rows = []
for ct, d in sig.groupby("cell_type"):
    for direction, genes in (("up", d[d.avg_log2FC > 0].gene), ("down", d[d.avg_log2FC < 0].gene)):
        genes = sorted(set(genes))
        if len(genes) < 10:
            continue
        for attempt in range(3):
            try:
                enr = gp.enrichr(gene_list=genes, gene_sets=libs, organism="human", outdir=None)
                break
            except Exception as e:  # network hiccup
                print("retry", ct, direction, e); time.sleep(5)
        else:
            continue
        r = enr.results
        r = r[r["Adjusted P-value"] < 0.05].copy()
        r["cell_type"], r["direction"], r["n_input_genes"] = ct, direction, len(genes)
        rows.append(r)
        print(f"{ct:30s} {direction:4s} genes={len(genes):4d} sig_terms={len(r)}")
        time.sleep(1)

res = pd.concat(rows)
res.to_csv("results/tables/enrichment_GO_KEGG_by_cell_type.csv", index=False)
top = (res.sort_values("Adjusted P-value").groupby(["cell_type", "direction", "Gene_set"]).head(5))
top[["cell_type", "direction", "Gene_set", "Term", "Overlap", "Adjusted P-value", "Genes"]] \
    .to_csv("results/tables/enrichment_top_terms.csv", index=False)
print(top[top.Gene_set.str.startswith("GO")][["cell_type", "direction", "Term", "Adjusted P-value"]].to_string())
