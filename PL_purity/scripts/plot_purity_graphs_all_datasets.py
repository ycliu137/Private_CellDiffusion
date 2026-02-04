"""
Plot neighbor purity for integration graphs built from X_fae and X_dif across datasets.
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

input_csvs = snakemake.input.csvs
output_csv = snakemake.output.csv
output_pdf = snakemake.output.pdf
dataset_names = snakemake.config.get("dataset_names") or {}

print(f"Reading {len(input_csvs)} CSV files")
frames = []
for csv_path in input_csvs:
    df = pd.read_csv(csv_path)
    if "dataset" not in df.columns:
        dataset_name = Path(csv_path).parent.name
        df["dataset"] = dataset_name
    frames.append(df)

if not frames:
    raise ValueError("No input CSV files provided")

combined = pd.concat(frames, ignore_index=True)
required_cols = ["dataset", "neighbor_purity_X_dif", "neighbor_purity_X_fae"]
missing_cols = [col for col in required_cols if col not in combined.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

# Sort by dataset name for consistent plotting
combined = combined.sort_values("dataset").reset_index(drop=True)

# Map dataset display names if provided
combined["dataset_display"] = combined["dataset"].map(lambda x: dataset_names.get(x, x))

# Save combined CSV
print(f"Saving combined CSV to: {output_csv}")
Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
combined.to_csv(output_csv, index=False)

# Plot grouped bar chart
print(f"Creating plot to: {output_pdf}")
fig, ax = plt.subplots(figsize=(12, 6))

x = np.arange(len(combined))
bar_width = 0.35

ax.bar(x - bar_width / 2, combined["neighbor_purity_X_dif"],
       width=bar_width, label="Before diffusion", color="#1f77b4")
ax.bar(x + bar_width / 2, combined["neighbor_purity_X_fae"],
       width=bar_width, label="After diffusion", color="#ff7f0e")

ax.set_xlabel("Datasets", fontsize=12, fontweight="bold")
ax.set_ylabel("Edge Accuracy", fontsize=12, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(combined["dataset_display"], rotation=45, ha="right")
ax.legend(fontsize=10, loc="best")
ax.grid(axis="y", alpha=0.3, linestyle="--")

plt.tight_layout()
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, format="pdf", dpi=300, bbox_inches="tight")
plt.close()

print("Plot saved successfully!")
