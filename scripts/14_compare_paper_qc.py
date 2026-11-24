import csv
from collections import defaultdict
from pathlib import Path

qc_file = Path("data/processed/GSE171213_HC_PD_qc_metrics.csv")
groups = defaultdict(lambda: {"before": 0, "after": 0})

with qc_file.open(newline="") as f:
    for row in csv.DictReader(f):
        donor = row["donor"]
        groups[donor]["before"] += 1

        # Cell filters stated in the original paper
        passes_paper_qc = (
            int(row["detected_genes"]) > 200
            and float(row["percent_mito"]) < 40
        )
        if passes_paper_qc:
            groups[donor]["after"] += 1

print("Donor  Condition       Before QC  After paper-style QC")
for donor in sorted(groups):
    condition = "Healthy" if donor.startswith("HC") else "Periodontitis"
    print(
        f"{donor:5}  {condition:14} "
        f"{groups[donor]['before']:9,}  {groups[donor]['after']:20,}"
    )

healthy_total = sum(
    v["after"] for d, v in groups.items() if d.startswith("HC")
)
pd_total = sum(
    v["after"] for d, v in groups.items() if d.startswith("PD")
)

print(f"\nHealthy total after these two filters: {healthy_total:,}")
print(f"Periodontitis total after these two filters: {pd_total:,}")
print("Paper reports: Healthy 14,552; Periodontitis 19,865")