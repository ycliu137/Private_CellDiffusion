"""
Plot SCIB evaluation results: bar plot comparing different graph building methods for aggregate scores
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import re

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

# Find all graph method rows (X_dif_{method})
print(f"\n=== Finding graph building methods in results ===")
method_rows = {}
for method_name in df.index:
    method_str = str(method_name)
    # Look for X_dif_{method} pattern
    if method_str.startswith('X_dif_'):
        method_key = method_str.replace('X_dif_', '')
        method_rows[method_key] = method_name
        print(f"Found method: {method_key} -> {method_name}")

if len(method_rows) == 0:
    raise ValueError("Could not find any X_dif methods in results table. Available methods: " + str(list(df.index)))

# Sort methods for consistent ordering
sorted_methods = sorted(method_rows.keys())
print(f"\nFound {len(sorted_methods)} graph building methods: {sorted_methods}")

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

# Extract scores for each method
scores_data = {}
for method_key in sorted_methods:
    method_name = method_rows[method_key]
    scores_data[method_key] = {}
    for score_name, col_name in found_scores.items():
        if col_name in df.columns:
            value = df.loc[method_name, col_name]
            if pd.isna(value):
                print(f"  Warning: {score_name} is NaN for {method_key}")
                value = 0.0
            scores_data[method_key][score_name] = float(value)
            print(f"  {method_key} - {score_name}: {value:.4f}")

# Create bar plot
print(f"\n=== Creating bar plot ===")

# Set up plot
fig, ax = plt.subplots(figsize=(14, 8))

# Set up bar positions
n_methods = len(sorted_methods)
n_scores = len(found_scores)
x = np.arange(n_methods)  # x positions for methods
width = 0.25  # Width of bars

# Colors for each aggregate score
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, Orange, Green
score_names = list(found_scores.keys())

# Plot bars for each aggregate score
for i, (score_name, col_name) in enumerate(found_scores.items()):
    values = [scores_data[method][score_name] for method in sorted_methods]
    offset = (i - 1) * width  # Center the bars
    ax.bar(x + offset, values, width, label=score_name, color=colors[i % len(colors)], 
           alpha=0.8, edgecolor='black', linewidth=1.2)

# Customize plot
ax.set_xlabel('Graph Building Method', fontsize=12, fontweight='bold')
ax.set_ylabel('Score Value', fontsize=12, fontweight='bold')
ax.set_title('SCIB Evaluation: Aggregate Scores Comparison - Graph Building Methods', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(sorted_methods, rotation=45, ha='right')  # Rotate labels for readability
ax.legend(loc='best', fontsize=11, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--', axis='y')

# Set y-axis limits dynamically based on max value
all_values = []
for method in sorted_methods:
    all_values.extend([scores_data[method][score] for score in found_scores.keys()])
min_val = min(all_values) if all_values else 0
max_val = max(all_values) if all_values else 1
ax.set_ylim(bottom=max(0, min_val - 0.1), top=min(1.0, max_val + 0.1))

# Add value labels on bars
for i, (score_name, col_name) in enumerate(found_scores.items()):
    values = [scores_data[method][score_name] for method in sorted_methods]
    offset = (i - 1) * width
    for j, (method, v) in enumerate(zip(sorted_methods, values)):
        ax.text(j + offset, v + 0.01, f'{v:.3f}', 
                ha='center', va='bottom', fontsize=8, fontweight='bold')

plt.tight_layout()

# Save figure
print(f"\n=== Saving plot ===")
print(f"Output file: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print(f"Bar plot saved successfully!")
print(f"  Compared methods: {sorted_methods}")
print(f"  Aggregate scores: {list(found_scores.keys())}")

print("\n=== Plotting complete! ===")

