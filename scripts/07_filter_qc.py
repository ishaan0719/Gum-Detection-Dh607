import csv
from pathlib import Path

import numpy as np
from scipy.sparse import load_npz, save_npz

processed_dir = Path("data/processed")

counts = load_npz(processed_dir / "GSE171213_HC_PD_counts.npz").tocsc()
genes = (processed_dir / "gene_names.txt").read_text().splitlines()

with open(processed_dir / "cell_names.txt") as handle:
    all_cell_names = handle.read().splitlines()

# These are the selected HC and PD cells, in matrix-column order.
with open(processed_dir / "GSE171213_HC_PD_cell_metadata.csv", newline="") as handle:
    metadata = list(csv.DictReader(handle))

with open(processed_dir / "GSE171213_HC_PD_qc_metrics.csv", newline="") as handle:
    qc_rows = list(csv.DictReader(handle))

if counts.shape[1] != len(metadata) or len(metadata) != len(qc_rows):
    raise ValueError("Counts, metadata, and QC rows do not match.")

# Published upper limits for each donor.
max_genes = {
    "HC1": 1500, "HC2": 1500, "HC3": 2000, "HC4": 2000,
    "PD1": 3000, "PD2": 2000, "PD3": 3000, "PD4": 3000, "PD5": 3000,
}
max_counts = {
    "HC1": 6000, "HC2": 6000, "HC3": 15000, "HC4": 10000,
    "PD1": 10000, "PD2": 10000, "PD3": 10000, "PD4": 15000, "PD5": 15000,
}

keep = np.zeros(counts.shape[1], dtype=bool)

for i, (meta, qc) in enumerate(zip(metadata, qc_rows)):
    donor = meta["donor"]
    n_genes = int(qc["detected_genes"])
    n_counts = int(qc["total_counts"])
    mito = float(qc["percent_mito"])

    keep[i] = (
        n_genes > 200
        and n_genes <= max_genes[donor]
        and n_counts > 3
        and n_counts <= max_counts[donor]
        and mito < 25
    )

filtered_counts = counts[:, keep].tocsr()
filtered_metadata = [row for row, retained in zip(metadata, keep) if retained]
filtered_qc = [row for row, retained in zip(qc_rows, keep) if retained]

save_npz(
    processed_dir / "GSE171213_HC_PD_filtered_counts.npz",
    filtered_counts,
)

with open(processed_dir / "GSE171213_HC_PD_filtered_metadata.csv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=filtered_metadata[0].keys())
    writer.writeheader()
    writer.writerows(filtered_metadata)

with open(processed_dir / "GSE171213_HC_PD_filtered_qc_metrics.csv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=filtered_qc[0].keys())
    writer.writeheader()
    writer.writerows(filtered_qc)

print(f"Cells before filtering: {len(keep):,}")
print(f"Cells retained: {int(keep.sum()):,}")
print(f"Cells removed: {int((~keep).sum()):,}")
print(f"Filtered matrix shape (genes × cells): {filtered_counts.shape}")

print("\nRetained cells per donor:")
for donor in sorted(max_genes):
    before = sum(row["donor"] == donor for row in metadata)
    after = sum(row["donor"] == donor for row in filtered_metadata)
    print(f"{donor}: {after:,} / {before:,}")