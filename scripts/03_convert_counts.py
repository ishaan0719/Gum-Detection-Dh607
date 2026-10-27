import csv
import gzip
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse import save_npz, vstack

project_dir = Path.cwd()
input_file = project_dir / "data" / "raw" / "GSE171213_AllSample.counts.tsv.gz"
output_dir = project_dir / "data" / "processed"
output_dir.mkdir(parents=True, exist_ok=True)

if not input_file.exists():
    raise FileNotFoundError(f"Count file not found: {input_file}")

print("Reading the count table and building a sparse matrix...")

cell_names = None
gene_names = []
matrix_chunks = []
chunk_rows = []
chunk_size = 200

with gzip.open(input_file, "rt", newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")

    # The first row contains "gene" followed by cell names.
    header = next(reader)
    if header[0] != "gene":
        raise ValueError(f"Expected first column to be 'gene'; found {header[0]!r}")

    cell_names = header[1:]
    number_of_cells = len(cell_names)

    for row_number, row in enumerate(reader, start=1):
        if len(row) != number_of_cells + 1:
            raise ValueError(
                f"Row {row_number + 1} has {len(row)} fields; "
                f"expected {number_of_cells + 1}."
            )

        gene_names.append(row[0])

        values = np.fromiter(
            (int(value) for value in row[1:]),
            dtype=np.int32,
            count=number_of_cells,
        )

        nonzero_columns = np.flatnonzero(values)
        nonzero_values = values[nonzero_columns]

        chunk_rows.append(
            sparse.csr_matrix(
                (
                    nonzero_values,
                    nonzero_columns,
                    np.array([0, len(nonzero_columns)]),
                ),
                shape=(1, number_of_cells),
                dtype=np.int32,
            )
        )

        if len(chunk_rows) == chunk_size:
            matrix_chunks.append(vstack(chunk_rows, format="csr"))
            chunk_rows = []

        if row_number % 2000 == 0:
            print(f"Read {row_number:,} genes...")

if chunk_rows:
    matrix_chunks.append(vstack(chunk_rows, format="csr"))

counts = vstack(matrix_chunks, format="csr")

save_npz(output_dir / "GSE171213_counts.npz", counts)

with open(output_dir / "gene_names.txt", "w") as handle:
    handle.write("\n".join(gene_names) + "\n")

with open(output_dir / "cell_names.txt", "w") as handle:
    handle.write("\n".join(cell_names) + "\n")

print("\nConversion complete.")
print(f"Genes: {counts.shape[0]:,}")
print(f"Cells: {counts.shape[1]:,}")
print(f"Nonzero entries: {counts.nnz:,}")
print(f"Saved files in: {output_dir}")