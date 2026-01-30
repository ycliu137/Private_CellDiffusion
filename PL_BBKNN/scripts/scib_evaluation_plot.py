"""
Bar plot of SCIB aggregate scores for BBKNN (single embedding X_bbknn).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

input_table = snakemake.input.table
output_pdf = snakemake.output.pdf

print("\n=== Loading SCIB results ===")
df = pd.read_csv(input_table, index_col=0)
print(f"Shape: {df.shape}, index: {list(df.index)}")

if "Metric Type" in df.index:
    df = df.drop("Metric Type")

# Expect a single embedding row (X_bbknn or similar)
embedding_rows = [i for i in df.index if str(i).startswith("X_")]
if not embedding_rows:
    embedding_rows = list(df.index)
row_name = embedding_rows[0]
display_name = "BBKNN" if "bbknn" in str(row_name).lower() else str(row_name)

agg = {}
for col in df.columns:
    c = str(col).lower()
    if col == "Total" or "total" in c:
        agg["Total"] = col
    elif "batch" in c and "correction" in c:
        agg["Batch correction"] = col
    elif "bio" in c and "conservation" in c:
        agg["Bio conservation"] = col
agg = {k: v for k, v in agg.items() if v is not None}
if not agg:
    raise ValueError("No aggregate score columns found.")

for c in agg.values():
    df[c] = pd.to_numeric(df[c], errors="coerce")

scores = [float(df.loc[row_name, agg[s]]) if not pd.isna(df.loc[row_name, agg[s]]) else 0.0 for s in agg]
labels = list(agg.keys())
colors = {"Total": "#4472C4", "Batch correction": "#70AD47", "Bio conservation": "#C55A11"}
col_list = [colors.get(l, "#808080") for l in labels]

fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(labels))
bars = ax.bar(x, scores, color=col_list, edgecolor="black", linewidth=1)
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Score")
ax.set_title(f"SCIB evaluation: {display_name}")
ax.set_ylim(0, 1.1)
for b, v in zip(bars, scores):
    ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.02, f"{v:.3f}", ha="center", va="bottom", fontsize=10)
plt.tight_layout()
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved: {output_pdf}")
