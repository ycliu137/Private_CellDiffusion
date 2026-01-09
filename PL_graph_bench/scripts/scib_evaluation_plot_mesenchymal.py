"""
Plot SCIB evaluation results for Mesenchymal lineage cells: bar plot comparing different graph building methods for aggregate scores
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import re

# Set publication-quality style parameters
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans', 'Liberation Sans']
rcParams['font.size'] = 11
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 14
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 10
rcParams['figure.titlesize'] = 14
rcParams['axes.linewidth'] = 1.0
rcParams['grid.linewidth'] = 0.5
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

# Get input and output from snakemake
input_table = snakemake.input.table
output_pdf = snakemake.output.pdf

print(f"\n=== Loading SCIB results (Mesenchymal) ===")
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

# Set up plot with dynamic size based on number of methods
n_methods = len(sorted_methods)
n_scores = len(found_scores)
fig_width = max(14, n_methods * 2.5)  # At least 2.5 units per method
fig, ax = plt.subplots(figsize=(fig_width, 8))

# Set up bar positions
x = np.arange(n_methods)  # x positions for methods
width = 0.28  # Width of bars (reduced for better spacing)

# Professional color palette
colors = {
    'Total': '#4472C4',  # Rich blue
    'Batch correction': '#70AD47',  # Green
    'Bio conservation': '#C55A11'  # Warm orange-red
}

score_names = list(found_scores.keys())

# Plot bars for each aggregate score
for i, (score_name, col_name) in enumerate(found_scores.items()):
    values = [scores_data[method][score_name] for method in sorted_methods]
    offset = (i - 1) * width  # Center the bars
    color = colors.get(score_name, plt.cm.tab10(i))
    ax.bar(x + offset, values, width, label=score_name, color=color, 
           alpha=0.9, edgecolor='white', linewidth=1.2, zorder=3)

# Customize plot
ax.set_xlabel('Graph Building Method', fontsize=12, fontweight='bold', labelpad=10)
ax.set_ylabel('Score', fontsize=12, fontweight='bold', labelpad=10)
ax.set_title('SCIB Evaluation (Mesenchymal Lineage): Aggregate Scores Comparison - Graph Building Methods', 
             fontsize=14, fontweight='bold', pad=15)
ax.set_xticks(x)
ax.set_xticklabels(sorted_methods, rotation=45, ha='right', fontsize=10)

# Improve grid styling
ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5, axis='y', zorder=0)
ax.set_axisbelow(True)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.0)
ax.spines['bottom'].set_linewidth(1.0)

# Legend
legend = ax.legend(loc='lower left', 
                   frameon=True, 
                   fancybox=False, 
                   shadow=False,
                   edgecolor='black',
                   facecolor='white',
                   framealpha=1.0,
                   borderpad=0.8,
                   handlelength=2.0,
                   handletextpad=0.5,
                   fontsize=10,
                   ncol=1)
legend.get_frame().set_linewidth(0.8)

# Set y-axis limits dynamically based on max value
all_values = []
for method in sorted_methods:
    all_values.extend([scores_data[method][score] for score in found_scores.keys()])
min_val = min(all_values) if all_values else 0
max_val = max(all_values) if all_values else 1
y_range = max_val - min_val
y_padding = y_range * 0.08 if y_range > 0 else 0.1
ax.set_ylim(bottom=max(0, min_val - y_padding), top=min(1.0, max_val + y_padding))

# Add minor ticks
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))

# Improve tick styling
ax.tick_params(axis='both', which='major', length=6, width=1.0)
ax.tick_params(axis='both', which='minor', length=3, width=0.5)

# Add value labels on bars (only if bars are tall enough)
for i, (score_name, col_name) in enumerate(found_scores.items()):
    values = [scores_data[method][score_name] for method in sorted_methods]
    offset = (i - 1) * width
    for j, (method, v) in enumerate(zip(sorted_methods, values)):
        # Only add label if bar is tall enough (>5% of max height)
        if v > (max_val - min_val) * 0.05:
            ax.text(j + offset, v + y_padding * 0.3, f'{v:.3f}', 
                    ha='center', va='bottom', fontsize=9, fontweight='bold',
                    zorder=4)

plt.tight_layout()

# Save figure
print(f"\n=== Saving plot ===")
print(f"Output file: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)

plt.savefig(output_pdf, 
           dpi=300, 
           bbox_inches='tight', 
           pad_inches=0.1,
           facecolor='white',
           edgecolor='none',
           format='pdf',
           transparent=False)

# Also save as PNG
png_output = str(output_pdf).replace('.pdf', '.png')
plt.savefig(png_output,
           dpi=300,
           bbox_inches='tight',
           pad_inches=0.1,
           facecolor='white',
           edgecolor='none',
           format='png',
           transparent=False)
print(f"Also saved PNG version: {png_output}")

plt.close()

print(f"Bar plot saved successfully!")
print(f"  Compared methods: {sorted_methods}")
print(f"  Aggregate scores: {list(found_scores.keys())}")

print("\n=== Plotting complete! ===")

