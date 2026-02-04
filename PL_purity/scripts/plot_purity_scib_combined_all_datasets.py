"""
Create a single-page PDF with two panels (side-by-side):
1) Neighbor purity bar plot (X_dif vs X_fae)
2) SCIB total score vs k_add line plot
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

# Inputs/outputs from snakemake
input_purity_csvs = snakemake.input.purity_csvs
input_scib_tables = snakemake.input.scib_tables
output_pdf = snakemake.output.pdf

dataset_names = snakemake.config.get("dataset_names") or {}

# Publication-style settings (font sizes only; text content unchanged)
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
})

# ----------------------
# Panel A: Neighbor purity bar plot
# ----------------------
frames = []
for csv_path in input_purity_csvs:
    df = pd.read_csv(csv_path)
    if "dataset" not in df.columns:
        dataset_name = Path(csv_path).parent.name
        df["dataset"] = dataset_name
    frames.append(df)

if not frames:
    raise ValueError("No purity CSV files provided")

purity_df = pd.concat(frames, ignore_index=True)
required_cols = ["dataset", "neighbor_purity_X_dif", "neighbor_purity_X_fae"]
missing_cols = [col for col in required_cols if col not in purity_df.columns]
if missing_cols:
    raise ValueError(f"Missing required columns in purity data: {missing_cols}")

purity_df = purity_df.sort_values("dataset").reset_index(drop=True)
purity_df["dataset_display"] = purity_df["dataset"].map(lambda x: dataset_names.get(x, x))

# ----------------------
# Panel B: SCIB total score line plot
# ----------------------
scib_series = []
for table_path in input_scib_tables:
    table_path = str(table_path)
    dataset_name = Path(table_path).parent.name
    display_name = dataset_names.get(dataset_name, dataset_name)

    df = pd.read_csv(table_path, index_col=0)
    if "Metric Type" in df.index:
        df = df.drop("Metric Type")

    total_col = None
    for col in df.columns:
        if col == "Total" or "total" in col.lower():
            total_col = col
            break
    if total_col is None:
        continue

    k_add_values = []
    for method_name in df.index:
        match = re.search(r"kadd(\d+)", str(method_name), re.IGNORECASE)
        k_add_values.append(int(match.group(1)) if match else None)

    df = df.copy()
    df["k_add"] = k_add_values
    df = df[df["k_add"].notna()].copy()
    if df.empty:
        continue

    df[total_col] = pd.to_numeric(df[total_col], errors='coerce')
    df = df.sort_values("k_add")

    scib_series.append((display_name, df["k_add"].values, df[total_col].values))

if not scib_series:
    raise ValueError("No datasets with valid Total score to plot")

# ----------------------
# Create combined figure
# ----------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)

# Panel A: bar chart
x = np.arange(len(purity_df))
bar_width = 0.35
ax1.bar(
    x - bar_width / 2,
    purity_df["neighbor_purity_X_dif"],
    width=bar_width,
    label="Before diffusion",
    color="#1f77b4",
)
ax1.bar(
    x + bar_width / 2,
    purity_df["neighbor_purity_X_fae"],
    width=bar_width,
    label="After diffusion",
    color="#ff7f0e",
)
ax1.set_xlabel("Datasets")
ax1.set_ylabel("Edge Accuracy")
ax1.set_xticks(x)
ax1.set_xticklabels(purity_df["dataset_display"], rotation=45, ha="right")
ax1.legend(loc="best")

# Panel B: line chart
for display_name, k_adds, totals in scib_series:
    ax2.plot(
        k_adds,
        totals,
        marker="o",
        linewidth=2,
        label=display_name,
    )
ax2.set_xlabel("Percentages of Added Random Edges")
ax2.set_ylabel("Total Score (scib-metrics)")
ax2.legend(loc="best", framealpha=0.9)

# Save output
print(f"Saving combined plot to: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches="tight")
plt.close()

print("Done.")
