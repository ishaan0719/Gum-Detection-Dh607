import csv
from collections import defaultdict
from pathlib import Path

qc_file = Path("data/processed/GSE171213_HC_PD_qc_metrics.csv")

by_donor = defaultdict(list)

with open(qc_file, newline="") as handle:
    for row in csv.DictReader(handle):
        by_donor[row["donor"]].append(row)

print("Per-donor QC summary")
print("Donor  Cells  Median genes  Median counts  Median mito%  % mito>25")

for donor in sorted(by_donor):
    rows = by_donor[donor]
    genes = [int(row["detected_genes"]) for row in rows]
    counts = [int(row["total_counts"]) for row in rows]
    mito = [float(row["percent_mito"]) for row in rows]

    print(
        f"{donor:5}  {len(rows):5,}  "
        f"{sorted(genes)[len(genes)//2]:12,}  "
        f"{sorted(counts)[len(counts)//2]:13,}  "
        f"{sorted(mito)[len(mito)//2]:12.1f}  "
        f"{100 * sum(value > 25 for value in mito) / len(mito):8.1f}%"
    )