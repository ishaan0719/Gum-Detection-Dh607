from pathlib import Path
import csv

import numpy as np
from scipy.sparse import load_npz

project_dir = Path.cwd()
processed_dir = project_dir / "data" / "processed"

counts = load_npz(processed_dir / "GSE171213_HC_PD_counts.npz").tocsc()
genes = (processed_dir / "gene_names.txt").read_text().splitlines()

metadata_file = processed_dir / "GSE171213_HC_PD_cell_metadata.csv"
with open(metadata_file, newline="") as handle:
    metadata = list(csv.DictReader(handle))

if counts.shape[0] != len(genes):
    raise ValueError("Gene-name count does not match the matrix rows.")
if counts.shape[1] != len(metadata):
    raise ValueError("Metadata rows do not match the matrix columns.")

# Metrics per cell
total_counts = np.asarray(counts.sum(axis=0)).ravel()
detected_genes = np.diff(counts.indptr)

# Human mitochondrial gene symbols usually start with MT-
mitochondrial_rows = [
    i for i, gene in enumerate(genes)
    if gene.upper().startswith("MT-")
]
mitochondrial_counts = np.asarray(
    counts[mitochondrial_rows, :].sum(axis=0)
).ravel()

percent_mito = np.divide(
    mitochondrial_counts * 100,
    total_counts,
    out=np.zeros_like(total_counts, dtype=float),
    where=total_counts > 0,
)

output_file = processed_dir / "GSE171213_HC_PD_qc_metrics.csv"
with open(output_file, "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "cell", "donor", "condition",
        "total_counts", "detected_genes", "percent_mito"
    ])
    for i, row in enumerate(metadata):
        writer.writerow([
            row["cell"],
            row["donor"],
            row["condition"],
            int(total_counts[i]),
            int(detected_genes[i]),
            round(float(percent_mito[i]), 4),
        ])

print(f"Saved QC metrics: {output_file}")
print(f"Mitochondrial genes found: {len(mitochondrial_rows)}")

for condition in ("Healthy", "Periodontitis"):
    indices = [
        i for i, row in enumerate(metadata)
        if row["condition"] == condition
    ]
    print(f"\n{condition} ({len(indices):,} cells)")
    for label, values in (
        ("Total counts", total_counts[indices]),
        ("Detected genes", detected_genes[indices]),
        ("Mitochondrial %", percent_mito[indices]),
    ):
        qs = np.percentile(values, [0, 25, 50, 75, 95, 100])
        print(
            f"{label}: min={qs[0]:.1f}, 25%={qs[1]:.1f}, "
            f"median={qs[2]:.1f}, 75%={qs[3]:.1f}, "
            f"95%={qs[4]:.1f}, max={qs[5]:.1f}"
        )