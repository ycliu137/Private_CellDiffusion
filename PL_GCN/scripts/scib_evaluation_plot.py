"""
Plot SCIB evaluation results: bar plot comparing CellDiffusion vs GCN aggregate scores
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np

# Get input and output from snakemake
input_table = snakemake.input.table
output_pdf = snakemake.output.pdf

print(f"\n=== Loading SCIB results ===")
print(f"Input file: {input_table}")

# Read CSV file
df = pd.read_csv(input_table, index_col=0)
print(f"Results table shape: {df.shape}")
print(f"Columns: {list(df.columns)}")
print(f"\nIndex (Embedding methods): {list(df.index)}")

# Check for "Metric Type" row and remove it if present
METRIC_TYPE_ROW = "Metric Type"
if METRIC_TYPE_ROW in df.index:
    print(f"\nRemoving '{METRIC_TYPE_ROW}' row from data")
    df = df.drop(METRIC_TYPE_ROW)

# Find CellDiffusion (X_dif) and GCN (X_gcn) rows
celldiffusion_row = None
gcn_row = None

print(f"\n=== Finding methods in results ===")
for method_name in df.index:
    method_str = str(method_name)
    # Look for X_dif but not X_dif_kadd* (to avoid matching stress test results)
    if 'X_dif' in method_str and 'X_dif_kadd' not in method_str:
        celldiffusion_row = method_name
        print(f"Found CellDiffusion: {method_name}")
    elif 'X_gcn' in method_str:
        gcn_row = method_name
        print(f"Found GCN: {method_name}")

if celldiffusion_row is None:
    # Try to find by exact match or partial match
    for idx in df.index:
        if 'X_dif' in str(idx) and 'kadd' not in str(idx):
            celldiffusion_row = idx
            print(f"Found CellDiffusion (alternative): {idx}")
            break
    if celldiffusion_row is None:
        raise ValueError("Could not find CellDiffusion (X_dif) in results table. Available methods: " + str(list(df.index)))

if gcn_row is None:
    # Try to find by exact match or partial match
    for idx in df.index:
        if 'X_gcn' in str(idx):
            gcn_row = idx
            print(f"Found GCN (alternative): {idx}")
            break
    if gcn_row is None:
        raise ValueError("Could not find GCN (X_gcn) in results table. Available methods: " + str(list(df.index)))

# Extract the three aggregate score columns
aggregate_score_cols = {
    'Total': None,
    'Batch correction': None,
    'Bio conservation': None
}

print(f"\n=== Extracting aggregate scores ===")
for col in df.columns:
    col_lower = col.lower()
    if col == 'Total' or 'total' in col_lower:
        aggregate_score_cols['Total'] = col
    elif 'batch' in col_lower and 'correction' in col_lower:
        aggregate_score_cols['Batch correction'] = col
    elif ('bio' in col_lower or 'conservation' in col_lower) and 'bio' in col_lower:
        aggregate_score_cols['Bio conservation'] = col

found_scores = {k: v for k, v in aggregate_score_cols.items() if v is not None}
print(f"Found aggregate score columns:")
for score_name, col_name in found_scores.items():
    print(f"  {score_name}: {col_name}")

if len(found_scores) == 0:
    raise ValueError("Could not find any aggregate score columns. Available columns: " + str(list(df.columns)))

# Convert numeric columns to float (in case they were read as strings)
for col_name in found_scores.values():
    df[col_name] = pd.to_numeric(df[col_name], errors='coerce')

# Extract scores for CellDiffusion and GCN
scores_data = {}
for method_name, method_label in [(celldiffusion_row, 'CellDiffusion'), (gcn_row, 'GCN')]:
    scores_data[method_label] = {}
    for score_name, col_name in found_scores.items():
        if col_name in df.columns:
            value = df.loc[method_name, col_name]
            if pd.isna(value):
                print(f"  Warning: {score_name} is NaN for {method_label}")
                value = 0.0
            scores_data[method_label][score_name] = float(value)
            print(f"  {method_label} - {score_name}: {value:.4f}")

# Create bar plot
print(f"\n=== Creating bar plot ===")

fig, ax = plt.subplots(figsize=(10, 6))

# Set up bar positions
x = np.arange(len(found_scores))
width = 0.35  # Width of bars

methods = list(scores_data.keys())
colors = ['#1f77b4', '#ff7f0e']  # Blue for CellDiffusion, Orange for GCN

# Plot bars
for i, method in enumerate(methods):
    values = [scores_data[method][score] for score in found_scores.keys()]
    offset = (i - 0.5) * width
    ax.bar(x + offset, values, width, label=method, color=colors[i], alpha=0.8, edgecolor='black', linewidth=1.2)

# Customize plot
ax.set_xlabel('Aggregate Score', fontsize=12, fontweight='bold')
ax.set_ylabel('Score Value', fontsize=12, fontweight='bold')
ax.set_title('SCIB Evaluation: CellDiffusion vs GCN - Aggregate Scores Comparison', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(list(found_scores.keys()), rotation=0, ha='center')
ax.legend(loc='best', fontsize=11, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--', axis='y')

# Set y-axis limits dynamically based on max value
all_values = []
for method in methods:
    all_values.extend([scores_data[method][score] for score in found_scores.keys()])
max_value = max(all_values) if all_values else 1.0
y_max = min(1.0, max_value * 1.1)  # Add 10% padding, cap at 1.0
ax.set_ylim(bottom=0, top=y_max)

# Add value labels on bars
for i, method in enumerate(methods):
    values = [scores_data[method][score] for score in found_scores.keys()]
    offset = (i - 0.5) * width
    for j, v in enumerate(values):
        ax.text(j + offset, v + 0.01, f'{v:.3f}', 
                ha='center', va='bottom', fontsize=9, fontweight='bold')

plt.tight_layout()

# Save figure
print(f"\n=== Saving plot ===")
print(f"Output file: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print(f"Bar plot saved successfully!")
print(f"  Compared methods: {methods}")
print(f"  Aggregate scores: {list(found_scores.keys())}")

print("\n=== Plotting complete! ===")

