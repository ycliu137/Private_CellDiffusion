"""
Plot SCIB evaluation results: aggregate scores vs k_add
Reads CSV file from scib_evaluation step and plots three aggregate scores as line plots.
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

# Check for "Metric Type" row and remove it if present
METRIC_TYPE_ROW = "Metric Type"
if METRIC_TYPE_ROW in df.index:
    print(f"\nRemoving '{METRIC_TYPE_ROW}' row from data")
    df = df.drop(METRIC_TYPE_ROW)

print(f"\nData after processing:")
print(df.head())

# Extract k_add values from Embedding column (index)
# Format: 'X_dif_kadd0', 'X_dif_kadd10', etc.
print(f"\n=== Extracting k_add values from method names ===")
k_add_values = []
for method_name in df.index:
    if 'kadd' in method_name.lower():
        match = re.search(r'kadd(\d+)', method_name, re.IGNORECASE)
        if match:
            k_add = int(match.group(1))
            k_add_values.append(k_add)
            print(f"  {method_name} -> k_add={k_add}")
        else:
            print(f"  Warning: Could not extract k_add from {method_name}")
            k_add_values.append(None)
    else:
        print(f"  Skipping {method_name} (no 'kadd' in name)")
        k_add_values.append(None)

# Add k_add column and filter valid rows
df['k_add'] = k_add_values
df = df[df['k_add'].notna()].copy()
df = df.sort_values('k_add')

print(f"\nFiltered data with {len(df)} rows")
print(f"\nData sorted by k_add:")
print(df[['k_add'] + [col for col in df.columns if col != 'k_add']].head(10))

# Extract the three aggregate score columns
# Expected columns: "Batch correction", "Bio conservation", "Total"
print(f"\n=== Extracting aggregate scores ===")
aggregate_score_cols = {
    'Total': None,
    'Batch correction': None,
    'Bio conservation': None
}

for col in df.columns:
    if col == 'k_add':
        continue
    col_lower = col.lower()
    if col == 'Total' or 'total' in col_lower:
        aggregate_score_cols['Total'] = col
    elif 'batch' in col_lower and 'correction' in col_lower:
        aggregate_score_cols['Batch correction'] = col
    elif ('bio' in col_lower or 'conservation' in col_lower) and 'bio' in col_lower:
        aggregate_score_cols['Bio conservation'] = col

# Check which scores were found
found_scores = {k: v for k, v in aggregate_score_cols.items() if v is not None}
print(f"Found aggregate score columns:")
for score_name, col_name in found_scores.items():
    print(f"  {score_name}: {col_name}")

if len(found_scores) == 0:
    raise ValueError("Could not find any aggregate score columns. Available columns: " + str(list(df.columns)))

# Plot aggregate scores vs k_add
print(f"\n=== Plotting aggregate scores vs k_add ===")

fig, ax = plt.subplots(figsize=(10, 6))

# Define plot styles for each score
plot_config = {
    'Total': {'color': '#1f77b4', 'marker': 'o', 'linestyle': '-', 'label': 'Total'},
    'Batch correction': {'color': '#ff7f0e', 'marker': 's', 'linestyle': '--', 'label': 'Batch correction'},
    'Bio conservation': {'color': '#2ca02c', 'marker': '^', 'linestyle': '-.', 'label': 'Bio conservation'}
}

# Collect all values to determine unified y-axis range
all_values = []

# Plot each aggregate score
for score_name, col_name in found_scores.items():
    if col_name in df.columns:
        # DataFrame is already sorted by k_add, so we can use values directly
        k_adds = df['k_add'].values
        values = df[col_name].values
        
        # Collect values for unified y-axis
        all_values.extend(values)
        
        # Debug: print the data being plotted
        print(f"\n  Plotting {score_name}:")
        print(f"    k_add values: {k_adds}")
        print(f"    {score_name} values: {values}")
        print(f"    {score_name} range: [{values.min():.4f}, {values.max():.4f}]")
        
        config = plot_config.get(score_name, {})
        ax.plot(k_adds, values,
                marker=config.get('marker', 'o'),
                linestyle=config.get('linestyle', '-'),
                color=config.get('color', '#1f77b4'),
                linewidth=2,
                markersize=8,
                label=config.get('label', score_name))

# Customize plot
ax.set_xlabel('k_add', fontsize=12, fontweight='bold')
ax.set_ylabel('Aggregate Score', fontsize=12, fontweight='bold')
ax.set_title('SCIB Evaluation: Aggregate Scores vs k_add', fontsize=14, fontweight='bold')
ax.legend(loc='best', fontsize=10, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--')

# Set x-axis limits with padding
if len(df) > 0:
    x_min, x_max = df['k_add'].min(), df['k_add'].max()
    ax.set_xlim(left=x_min - 5, right=x_max + 5)

# Set unified y-axis range for all aggregate scores
if all_values:
    y_min = min(all_values)
    y_max = max(all_values)
    # Add padding (5% on each side)
    y_range = y_max - y_min
    y_padding = y_range * 0.05
    ax.set_ylim(bottom=y_min - y_padding, top=y_max + y_padding)
    print(f"\n  Unified y-axis range: [{y_min - y_padding:.4f}, {y_max + y_padding:.4f}]")

plt.tight_layout()

# Save figure
print(f"\n=== Saving plot ===")
print(f"Output file: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved successfully!")
print(f"  Plotted {len(found_scores)} aggregate scores:")
for score_name in found_scores.keys():
    print(f"    - {score_name}")

print("\n=== Plotting complete! ===")
