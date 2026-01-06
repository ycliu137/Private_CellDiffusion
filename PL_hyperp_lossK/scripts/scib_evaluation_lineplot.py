"""
Plot SCIB evaluation results: line plot comparing multiple datasets across loss_k values
Shows three aggregate scores (Total, Batch correction, Bio conservation) as lines for each dataset
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
rcParams['lines.linewidth'] = 2.0
rcParams['patch.linewidth'] = 0.5
rcParams['xtick.major.width'] = 1.0
rcParams['ytick.major.width'] = 1.0
rcParams['xtick.minor.width'] = 0.5
rcParams['ytick.minor.width'] = 0.5
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

# Get input and output from snakemake
input_tables = snakemake.input.tables  # List of CSV files
output_pdf = snakemake.output.pdf
dataset_dict = snakemake.params.dataset_dict  # Dictionary mapping dataset keys to display names

print(f"\n=== Loading SCIB results from {len(input_tables)} datasets ===")

# Store data for all datasets
all_datasets_data = {}

# Process each dataset
for table_file in input_tables:
    print(f"\nProcessing: {table_file}")
    
    # Extract dataset key from file path
    # Path format: .../output_dir/{dataset}/scib_results_table.csv
    path_parts = Path(table_file).parts
    dataset_key = None
    for part in reversed(path_parts):
        if part.endswith('.csv'):
            continue
        # Find the dataset directory name
        dataset_key = part
        break
    
    if dataset_key is None:
        print(f"  Warning: Could not extract dataset key from path {table_file}, skipping")
        continue
    
    # Get display name from dictionary, or use key as fallback
    display_name = dataset_dict.get(dataset_key, dataset_key)
    print(f"  Dataset key: {dataset_key}, Display name: {display_name}")
    
    # Read CSV file
    try:
        df = pd.read_csv(table_file, index_col=0)
        print(f"  Results table shape: {df.shape}")
    except Exception as e:
        print(f"  Error reading {table_file}: {e}, skipping")
        continue
    
    # Check for "Metric Type" row and remove it if present
    METRIC_TYPE_ROW = "Metric Type"
    if METRIC_TYPE_ROW in df.index:
        df = df.drop(METRIC_TYPE_ROW)
    
    # Find all loss_k rows (X_dif_lossk{loss_k})
    loss_k_methods = {}
    for method_name in df.index:
        method_str = str(method_name)
        match = re.search(r'X_dif_lossk(\d+)', method_str)
        if match:
            loss_k = int(match.group(1))
            loss_k_methods[loss_k] = method_name
    
    if len(loss_k_methods) == 0:
        print(f"  Warning: No X_dif_lossk methods found in {table_file}, skipping")
        continue
    
    # Sort by loss_k value
    sorted_loss_k = sorted(loss_k_methods.keys())
    print(f"  Found {len(sorted_loss_k)} loss_k values: {sorted_loss_k}")
    
    # Extract the three aggregate score columns
    aggregate_score_cols = {
        'Total': None,
        'Batch correction': None,
        'Bio conservation': None
    }
    
    for col in df.columns:
        col_lower = col.lower()
        if col == 'Total' or 'total' in col_lower:
            aggregate_score_cols['Total'] = col
        elif 'batch' in col_lower and 'correction' in col_lower:
            aggregate_score_cols['Batch correction'] = col
        elif ('bio' in col_lower or 'conservation' in col_lower) and 'bio' in col_lower:
            aggregate_score_cols['Bio conservation'] = col
    
    found_scores = {k: v for k, v in aggregate_score_cols.items() if v is not None}
    if len(found_scores) == 0:
        print(f"  Warning: No aggregate score columns found in {table_file}, skipping")
        continue
    
    # Convert numeric columns to float
    for col_name in found_scores.values():
        df[col_name] = pd.to_numeric(df[col_name], errors='coerce')
    
    # Extract scores for each loss_k value
    scores_data = {}
    for loss_k in sorted_loss_k:
        method_name = loss_k_methods[loss_k]
        scores_data[loss_k] = {}
        for score_name, col_name in found_scores.items():
            if col_name in df.columns:
                value = df.loc[method_name, col_name]
                if pd.isna(value):
                    print(f"    Warning: {score_name} is NaN for loss_k={loss_k}")
                    value = 0.0
                scores_data[loss_k][score_name] = float(value)
    
    # Store data for this dataset
    all_datasets_data[display_name] = {
        'loss_k_values': sorted_loss_k,
        'scores_data': scores_data,
        'found_scores': found_scores
    }
    print(f"  Successfully loaded data for {display_name}")

if len(all_datasets_data) == 0:
    raise ValueError("No valid datasets found to plot")

print(f"\n=== Creating multi-dataset line plot ===")
print(f"Total datasets to plot: {len(all_datasets_data)}")

# Get all unique loss_k values across all datasets
all_loss_k_values = set()
for dataset_name, data in all_datasets_data.items():
    all_loss_k_values.update(data['loss_k_values'])

sorted_all_loss_k = sorted(all_loss_k_values)
print(f"All loss_k values: {sorted_all_loss_k}")

# Get score names (should be the same for all datasets)
score_names = None
for dataset_name, data in all_datasets_data.items():
    if score_names is None:
        score_names = list(data['found_scores'].keys())
        break

if score_names is None:
    raise ValueError("No score names found")

print(f"Score types: {score_names}")

# Create figure with subplots for each score type
n_scores = len(score_names)
fig, axes = plt.subplots(1, n_scores, figsize=(6 * n_scores, 6))

if n_scores == 1:
    axes = [axes]
else:
    axes = axes.flatten()

# Color palette for datasets - use distinct colors
dataset_colors = plt.cm.tab10(np.linspace(0, 1, len(all_datasets_data)))
# Alternative: use specific colors
color_list = ['#4472C4', '#70AD47', '#C55A11', '#7030A0', '#FFC000', '#ED7D31', '#5B9BD5', '#A5A5A5']

# Color and marker styles for each dataset
for idx, dataset_name in enumerate(sorted(all_datasets_data.keys())):
    color = color_list[idx % len(color_list)]
    marker = ['o', 's', '^', 'D', 'v', 'p', '*', 'h'][idx % 8]
    linestyle = '-'
    
    data = all_datasets_data[dataset_name]
    loss_k_values = data['loss_k_values']
    scores_data = data['scores_data']
    
    # Plot each score type in its own subplot
    for score_idx, score_name in enumerate(score_names):
        ax = axes[score_idx]
        
        # Get values for this score
        values = []
        valid_loss_k = []
        for loss_k in sorted_all_loss_k:
            if loss_k in scores_data and score_name in scores_data[loss_k]:
                values.append(scores_data[loss_k][score_name])
                valid_loss_k.append(loss_k)
        
        if len(values) > 0:
            ax.plot(valid_loss_k, values, 
                   marker=marker, 
                   linestyle=linestyle,
                   color=color,
                   label=dataset_name,
                   linewidth=2.0,
                   markersize=6,
                   markeredgewidth=1.0,
                   markeredgecolor='white',
                   alpha=0.9)
        
        # Customize subplot
        ax.set_xlabel('loss_k', fontsize=12, fontweight='bold', labelpad=8)
        ax.set_ylabel(score_name, fontsize=12, fontweight='bold', labelpad=8)
        ax.set_title(score_name, fontsize=14, fontweight='bold', pad=10)
        
        # Improve grid styling
        ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5, axis='both', zorder=0)
        ax.set_axisbelow(True)
        
        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.0)
        ax.spines['bottom'].set_linewidth(1.0)
        
        # Set x-axis limits with padding
        if len(valid_loss_k) > 0:
            x_min = min(valid_loss_k)
            x_max = max(valid_loss_k)
            x_padding = (x_max - x_min) * 0.05 if (x_max - x_min) > 0 else 1
            ax.set_xlim(left=max(0, x_min - x_padding), right=x_max + x_padding)
        
        # Set y-axis limits
        all_values_this_score = []
        for dataset_name2, data2 in all_datasets_data.items():
            for loss_k in data2['loss_k_values']:
                if loss_k in data2['scores_data'] and score_name in data2['scores_data'][loss_k]:
                    all_values_this_score.append(data2['scores_data'][loss_k][score_name])
        
        if len(all_values_this_score) > 0:
            min_val = min(all_values_this_score)
            max_val = max(all_values_this_score)
            y_padding = (max_val - min_val) * 0.08 if (max_val - min_val) > 0 else 0.1
            ax.set_ylim(bottom=max(0, min_val - y_padding), 
                       top=min(1.0, max_val + y_padding))
        
        # Add minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))
        
        # Improve tick styling
        ax.tick_params(axis='both', which='major', length=6, width=1.0)
        ax.tick_params(axis='both', which='minor', length=3, width=0.5)
        ax.tick_params(axis='x', pad=5)
        ax.tick_params(axis='y', pad=5)
        
        # Add legend to first subplot only
        if score_idx == 0:
            legend = ax.legend(loc='best',
                             frameon=True,
                             fancybox=False,
                             shadow=False,
                             edgecolor='black',
                             facecolor='white',
                             framealpha=1.0,
                             borderpad=0.8,
                             handlelength=2.0,
                             handletextpad=0.5,
                             columnspacing=1.0,
                             fontsize=10,
                             ncol=1)
            legend.get_frame().set_linewidth(0.8)

# Tight layout
plt.tight_layout(pad=2.0)

# Save figure
print(f"\n=== Saving line plot ===")
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

print(f"\nLine plot saved successfully!")
print(f"  Datasets compared: {len(all_datasets_data)}")
print(f"  Score types: {score_names}")
print("\n=== Plotting complete! ===")

