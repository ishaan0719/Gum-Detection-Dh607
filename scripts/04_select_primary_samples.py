from pathlib import Path
import csv

import numpy as np
from scipy.sparse import load_npz, save_npz

project_dir = Path.cwd()
processed_dir = project_dir / "data" / "processed"

counts = load_npz(processed_dir / "GSE171213_counts.npz")
cell_names = (processed_dir / "cell_names.txt").read_text().splitlines()

if counts.shape[1] != len(cell_names):
    raise ValueError("Cell-name count does not match the matrix columns.")

# Each cell name starts with its donor ID, e.g. HC1_C1 or PD3_C25.
donors_to_keep = {f"HC{i}" for i in range(1, 5)}
donors_to_keep.update({f"PD{i}" for i in range(1, 6)})

donors = [name.split("_", 1)[0] for name in cell_names]
keep_indices = [
    i for i, donor in enumerate(donors)
    if donor in donors_to_keep
]

if not keep_indices:
    raise ValueError("No HC1–HC4 or PD1–PD5 cells matched the names.")

selected_counts = counts[:, keep_indices].tocsr()
selected_cells = [cell_names[i] for i in keep_indices]

output_matrix = processed_dir / "GSE171213_HC_PD_counts.npz"
save_npz(output_matrix, selected_counts)

metadata_file = processed_dir / "GSE171213_HC_PD_cell_metadata.csv"
with open(metadata_file, "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["cell", "donor", "condition"])
    for cell in selected_cells:
        donor = cell.split("_", 1)[0]
        condition = "Healthy" if donor.startswith("HC") else "Periodontitis"
        writer.writerow([cell, donor, condition])

print(f"Selected cells: {len(selected_cells):,}")
print(f"Matrix shape (genes × cells): {selected_counts.shape}")
print(f"Saved: {output_matrix}")
print(f"Saved: {metadata_file}")

print("\nCells per donor:")
for donor in sorted(donors_to_keep):
    print(f"{donor}: {sum(d == donor for d in (selected_cells_cell.split('_', 1)[0] for selected_cells_cell in selected_cells)):,}")