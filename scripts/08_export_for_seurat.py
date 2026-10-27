from pathlib import Path
import gzip

from scipy.sparse import load_npz
from scipy.io import mmwrite

processed_dir = Path("data/processed")
matrix = load_npz(
    processed_dir / "GSE171213_HC_PD_filtered_counts.npz"
).tocsr()

output_dir = processed_dir / "seurat_input"
output_dir.mkdir(parents=True, exist_ok=True)

# Matrix Market stores this as genes × cells, the orientation Seurat expects.
with gzip.open(output_dir / "counts.mtx.gz", "wb") as handle:
    mmwrite(handle, matrix)

metadata = processed_dir / "GSE171213_HC_PD_filtered_metadata.csv"
(output_dir / "metadata.csv").write_bytes(metadata.read_bytes())

genes = (processed_dir / "gene_names.txt").read_text().splitlines()
cells = [
    row.strip()
    for row in (processed_dir / "GSE171213_HC_PD_filtered_metadata.csv")
        .read_text()
        .splitlines()[1:]
]

# Save one gene and cell name per line, in matrix order.
(output_dir / "genes.txt").write_text("\n".join(genes) + "\n")
(output_dir / "cells.txt").write_text(
    "\n".join(row.split(",")[0] for row in cells) + "\n"
)

print(f"Matrix shape (genes × cells): {matrix.shape}")
print(f"Exported Seurat input files to: {output_dir}")