"""
Plot SCIB total score vs k_add for all datasets on the same figure.
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import re

input_tables = snakemake.input.tables
output_pdf = snakemake.output.pdf
dataset_names = snakemake.config.get("dataset_names") or {}

print(f"Reading {len(input_tables)} SCIB result tables")

fig, ax = plt.subplots(figsize=(10, 6))
plotted = 0

for table_path in input_tables:
    table_path = str(table_path)
    dataset_name = Path(table_path).parent.name
    display_name = dataset_names.get(dataset_name, dataset_name)
    df = pd.read_csv(table_path, index_col=0)

    # Remove metric type row if present
    if "Metric Type" in df.index:
        df = df.drop("Metric Type")

    # Identify total column
    total_col = None
    for col in df.columns:
        if col == "Total" or "total" in col.lower():
            total_col = col
            break
    if total_col is None:
        print(f"Warning: No Total column in {table_path}, skipping")
        continue

    # Extract k_add from index (e.g., X_dif_kadd10)
    k_add_values = []
    for method_name in df.index:
        match = re.search(r"kadd(\d+)", str(method_name), re.IGNORECASE)
        k_add_values.append(int(match.group(1)) if match else None)

    df = df.copy()
    df["k_add"] = k_add_values
    df = df[df["k_add"].notna()].copy()
    if df.empty:
        print(f"Warning: No k_add values found in {table_path}, skipping")
        continue

    df[total_col] = pd.to_numeric(df[total_col], errors='coerce')
    df = df.sort_values("k_add")

    ax.plot(
        df["k_add"].values,
        df[total_col].values,
        marker='o',
        linewidth=2,
        label=display_name
    )
    plotted += 1

if plotted == 0:
    raise ValueError("No datasets with valid Total score to plot")

ax.set_xlabel("Percentages of Added Random Edges", fontsize=12, fontweight='bold')
ax.set_ylabel("Total Score (scib-metrics)", fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3, linestyle='--')
ax.legend(loc='best', fontsize=8, framealpha=0.9)

plt.tight_layout()
print(f"Saving plot to: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print("Done.")
