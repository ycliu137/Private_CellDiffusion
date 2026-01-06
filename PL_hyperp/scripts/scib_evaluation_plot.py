"""
Plot SCIB evaluation results: bar plot comparing different k_mnn values for aggregate scores
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
# Use Arial/Helvetica font family for better compatibility
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
rcParams['lines.linewidth'] = 1.5
rcParams['patch.linewidth'] = 0.5
rcParams['xtick.major.width'] = 1.0
rcParams['ytick.major.width'] = 1.0
rcParams['xtick.minor.width'] = 0.5
rcParams['ytick.minor.width'] = 0.5
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

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

# Find all k_mnn rows (X_dif_kmnn{k_mnn})
print(f"\n=== Finding k_mnn methods in results ===")
k_mnn_methods = {}
for method_name in df.index:
    method_str = str(method_name)
    # Look for X_dif_kmnn{k_mnn} pattern
    match = re.search(r'X_dif_kmnn(\d+)', method_str)
    if match:
        k_mnn = int(match.group(1))
        k_mnn_methods[k_mnn] = method_name
        print(f"Found k_mnn={k_mnn}: {method_name}")

if len(k_mnn_methods) == 0:
    raise ValueError("Could not find any X_dif_kmnn methods in results table. Available methods: " + str(list(df.index)))

# Sort by k_mnn value
sorted_k_mnn = sorted(k_mnn_methods.keys())
print(f"\nFound {len(sorted_k_mnn)} k_mnn values: {sorted_k_mnn}")

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

# Extract scores for each k_mnn value
scores_data = {}
for k_mnn in sorted_k_mnn:
    method_name = k_mnn_methods[k_mnn]
    scores_data[k_mnn] = {}
    for score_name, col_name in found_scores.items():
        if col_name in df.columns:
            value = df.loc[method_name, col_name]
            if pd.isna(value):
                print(f"  Warning: {score_name} is NaN for k_mnn={k_mnn}")
                value = 0.0
            scores_data[k_mnn][score_name] = float(value)
            print(f"  k_mnn={k_mnn} - {score_name}: {value:.4f}")

# Create bar plot
print(f"\n=== Creating bar plot ===")

# Calculate figure size based on number of k_mnn values
# Use a wider aspect ratio suitable for publication
n_k_mnn = len(sorted_k_mnn)
n_scores = len(found_scores)
base_width = max(12, n_k_mnn * 1.2)  # Dynamic width based on number of bars
fig_height = 6  # Fixed height for consistent appearance
fig, ax = plt.subplots(figsize=(base_width, fig_height))

# Set up bar positions with better spacing
x = np.arange(n_k_mnn)  # x positions for k_mnn values
bar_width = 0.28  # Slightly wider bars for better visibility
bar_spacing = 0.02  # Small gap between groups

# Professional color palette - alternative scheme for publication
# Elegant and colorblind-friendly palette with good contrast
colors_palette = {
    'Total': '#4472C4',           # Rich blue
    'Batch correction': '#70AD47', # Forest green
    'Bio conservation': '#C55A11'  # Warm orange-red
}
# Fallback colors if score names don't match
fallback_colors = ['#4472C4', '#70AD47', '#C55A11', '#7030A0', '#FFC000']
score_names = list(found_scores.keys())

# Plot bars for each aggregate score with improved styling
bar_handles = []
for i, (score_name, col_name) in enumerate(found_scores.items()):
    values = [scores_data[k_mnn][score_name] for k_mnn in sorted_k_mnn]
    offset = (i - (n_scores - 1) / 2) * (bar_width + bar_spacing)  # Center the bars
    
    # Get color for this score
    color = colors_palette.get(score_name, fallback_colors[i % len(fallback_colors)])
    
    bars = ax.bar(x + offset, values, bar_width, 
                  label=score_name, 
                  color=color, 
                  alpha=0.9,  # Slightly more opaque
                  edgecolor='white',  # White edges for better separation
                  linewidth=1.0,  # Thinner edge for cleaner look
                  capsize=0)  # No caps on error bars
    
    bar_handles.append(bars)
    
    # Add value labels on bars (only if bar is tall enough)
    for j, (k_mnn, v) in enumerate(zip(sorted_k_mnn, values)):
        # Only label if the bar is tall enough to fit text
        y_pos = v + 0.005  # Small offset above bar
        if v > (max(values) * 0.05):  # Only label if bar is >5% of max height
            ax.text(j + offset, y_pos, f'{v:.3f}', 
                   ha='center', va='bottom', 
                   fontsize=9, 
                   fontweight='normal',
                   color='black')

# Customize plot with publication-quality styling
ax.set_xlabel('k_mnn', fontsize=12, fontweight='bold', labelpad=8)
ax.set_ylabel('Score', fontsize=12, fontweight='bold', labelpad=8)
# Remove title for cleaner look (can add back if needed)
# ax.set_title('SCIB Evaluation: Aggregate Scores vs k_mnn', fontsize=14, fontweight='bold', pad=15)

# Set x-axis
ax.set_xticks(x)
ax.set_xticklabels(sorted_k_mnn, fontsize=10, fontweight='normal')

# Improve legend styling and placement
legend = ax.legend(loc='lower left', 
                   frameon=True, 
                   fancybox=False, 
                   shadow=False,
                   edgecolor='black',
                   facecolor='white',
                   framealpha=1.0,
                   borderpad=0.8,
                   handlelength=1.5,
                   handletextpad=0.5,
                   columnspacing=1.0,
                   fontsize=10,
                   ncol=1)  # Single column for cleaner look

# Set legend frame properties
legend.get_frame().set_linewidth(0.8)

# Improve grid styling - subtle and professional
ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5, axis='y', zorder=0)
ax.set_axisbelow(True)  # Grid behind bars

# Remove top and right spines for cleaner look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.0)
ax.spines['bottom'].set_linewidth(1.0)

# Set y-axis limits with better padding
all_values = []
for k_mnn in sorted_k_mnn:
    all_values.extend([scores_data[k_mnn][score] for score in found_scores.keys()])
min_val = min(all_values) if all_values else 0
max_val = max(all_values) if all_values else 1
y_padding = (max_val - min_val) * 0.08  # 8% padding
ax.set_ylim(bottom=max(0, min_val - y_padding), 
            top=min(1.0, max_val + y_padding))

# Ensure y-axis starts at 0 for better comparison (if all values are positive)
if min_val >= 0:
    ax.set_ylim(bottom=0, top=min(1.0, max_val + y_padding))

# Add minor ticks for better readability
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))

# Improve tick styling
ax.tick_params(axis='both', which='major', length=6, width=1.0)
ax.tick_params(axis='both', which='minor', length=3, width=0.5)
ax.tick_params(axis='x', pad=5)
ax.tick_params(axis='y', pad=5)

# Tight layout with proper padding
plt.tight_layout(pad=2.0)

# Save figure with publication-quality settings
print(f"\n=== Saving plot ===")
print(f"Output file: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)

# Save with high resolution and proper formatting
plt.savefig(output_pdf, 
           dpi=300, 
           bbox_inches='tight', 
           pad_inches=0.1,
           facecolor='white',
           edgecolor='none',
           format='pdf',
           transparent=False)

# Also save as PNG for preview (optional)
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
print(f"  Compared k_mnn values: {sorted_k_mnn}")
print(f"  Aggregate scores: {list(found_scores.keys())}")

print("\n=== Plotting complete! ===")

