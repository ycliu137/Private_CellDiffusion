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

# Find all graph method rows (X_dif_{method} and X_gcn_{method})
print(f"\n=== Finding integration methods in results ===")
method_rows = {}  # {(method, integration_type): row_name}
for method_name in df.index:
    method_str = str(method_name)
    # Look for X_dif_{method} pattern (CellDiffusion)
    if method_str.startswith('X_dif_'):
        method_key = method_str.replace('X_dif_', '')
        method_rows[(method_key, 'CellDiffusion')] = method_name
        print(f"Found CellDiffusion method: {method_key} -> {method_name}")
    # Look for X_gcn_{method} pattern (GCN)
    elif method_str.startswith('X_gcn_'):
        method_key = method_str.replace('X_gcn_', '')
        method_rows[(method_key, 'GCN')] = method_name
        print(f"Found GCN method: {method_key} -> {method_name}")

if len(method_rows) == 0:
    raise ValueError("Could not find any X_dif or X_gcn methods in results table. Available methods: " + str(list(df.index)))

# Extract unique graph building methods and integration types
all_method_keys = sorted(set(method for method, _ in method_rows.keys()))
integration_types = ['CellDiffusion', 'GCN']
print(f"\nFound {len(method_rows)} integration methods across {len(all_method_keys)} graph building methods")
print(f"Graph building methods: {all_method_keys}")
print(f"Integration types: {integration_types}")

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

# Extract scores for each method and integration type
scores_data = {}  # {(method, integration_type): {score_name: value}}
for (method_key, integration_type), method_name in method_rows.items():
    key = (method_key, integration_type)
    scores_data[key] = {}
    for score_name, col_name in found_scores.items():
        if col_name in df.columns:
            value = df.loc[method_name, col_name]
            if pd.isna(value):
                print(f"  Warning: {score_name} is NaN for {method_key} ({integration_type})")
                value = 0.0
            scores_data[key][score_name] = float(value)
            print(f"  {method_key} ({integration_type}) - {score_name}: {value:.4f}")

# Create bar plot
print(f"\n=== Creating bar plot ===")

# Set up plot
fig, ax = plt.subplots(figsize=(16, 8))

# Set up bar positions
# For each graph method, we'll have bars for CellDiffusion and GCN, grouped by score type
n_methods = len(all_method_keys)
n_integration_types = len(integration_types)
n_scores = len(found_scores)
width = 0.12  # Width of bars (smaller to fit more bars)

# Colors for each aggregate score
score_colors = {'Total': '#4472C4', 'Batch correction': '#70AD47', 'Bio conservation': '#C55A11'}
# Different patterns/shades for integration types
integration_alpha = {'CellDiffusion': 0.9, 'GCN': 0.7}
integration_edge = {'CellDiffusion': 'black', 'GCN': 'gray'}

x = np.arange(n_methods)  # x positions for methods
score_names = list(found_scores.keys())

# Plot bars: for each method, show CellDiffusion and GCN side by side for each score
bar_idx = 0

for i, score_name in enumerate(score_names):
    for j, integration_type in enumerate(integration_types):
        values = []
        for method_key in all_method_keys:
            key = (method_key, integration_type)
            if key in scores_data:
                values.append(scores_data[key][score_name])
            else:
                values.append(0.0)
        
        offset = bar_idx * width - (n_scores * n_integration_types * width) / 2 + width / 2
        color = score_colors.get(score_name, '#808080')
        
        ax.bar(x + offset, values, width, 
               alpha=integration_alpha[integration_type],
               edgecolor=integration_edge[integration_type],
               linewidth=1.0,
               color=color,
               label=f'{score_name} ({integration_type})')
        
        bar_idx += 1

# Customize plot
ax.set_xlabel('Graph Building Method', fontsize=12, fontweight='bold')
ax.set_ylabel('Score Value', fontsize=12, fontweight='bold')
ax.set_title('SCIB Evaluation: Aggregate Scores Comparison (CellDiffusion & GCN)', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(all_method_keys, rotation=45, ha='right')
ax.legend(loc='best', fontsize=10, framealpha=0.9, ncol=3)
ax.grid(True, alpha=0.3, linestyle='--', axis='y')

# Set y-axis limits dynamically based on max value
all_values = []
for key in scores_data.keys():
    all_values.extend(list(scores_data[key].values()))
min_val = min(all_values) if all_values else 0
max_val = max(all_values) if all_values else 1
ax.set_ylim(bottom=max(0, min_val - 0.1), top=min(1.0, max_val + 0.1))

# Add value labels on bars (only if bars are tall enough)
for i, score_name in enumerate(score_names):
    for j, integration_type in enumerate(integration_types):
        values = []
        for method_key in all_method_keys:
            key = (method_key, integration_type)
            values.append(scores_data.get(key, {}).get(score_name, 0.0))
        
        offset = (i * n_integration_types + j) * width - (n_scores * n_integration_types * width) / 2 + width / 2
        for k, (method_key, v) in enumerate(zip(all_method_keys, values)):
            if v > 0.05:  # Only label if bar is tall enough
                ax.text(k + offset, v + 0.01, f'{v:.3f}', 
                        ha='center', va='bottom', fontsize=7, fontweight='bold')

plt.tight_layout()

# Save figure
print(f"\n=== Saving plot ===")
print(f"Output file: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print(f"Bar plot saved successfully!")
print(f"  Compared graph methods: {all_method_keys}")
print(f"  Integration types: {integration_types}")
print(f"  Aggregate scores: {list(found_scores.keys())}")

print("\n=== Plotting complete! ===")

