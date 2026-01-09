"""
Plot SCIB evaluation results for Mesenchymal lineage cells: line plot comparing CellDiffusion vs GCN aggregate scores across network layers
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
rcParams['lines.linewidth'] = 2.0
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

# Get input/output from snakemake
input_csv = snakemake.input.csv
output_pdf = snakemake.output.pdf

print(f"\n=== Loading SCIB results (Mesenchymal) ===")
print(f"Input file: {input_csv}")

df = pd.read_csv(input_csv, index_col=0)

print(f"\n=== Data summary ===")
print(f"Results table shape: {df.shape}")
print(f"\nColumn names: {list(df.columns)}")
print(f"\nIndex (Embedding names): {list(df.index)}")

# Extract network layers and method from embedding names
# Format: X_dif_nsteps{N} or X_gcn_nlayers{N}
network_layers_list = []
methods_list = []
valid_indices = []

for idx in df.index:
    # Try to match CellDiffusion pattern: X_dif_nsteps{N}
    match_dif = re.search(r'X_dif_nsteps(\d+)', str(idx))
    # Try to match GCN pattern: X_gcn_nlayers{N}
    match_gcn = re.search(r'X_gcn_nlayers(\d+)', str(idx))
    
    if match_dif:
        network_layers = int(match_dif.group(1))
        method = 'CellDiffusion'
        network_layers_list.append(network_layers)
        methods_list.append(method)
        valid_indices.append(idx)
    elif match_gcn:
        network_layers = int(match_gcn.group(1))
        method = 'GCN'
        network_layers_list.append(network_layers)
        methods_list.append(method)
        valid_indices.append(idx)
    else:
        print(f"Warning: Could not parse embedding name: {idx}, skipping")

# Filter dataframe to only valid indices
if len(valid_indices) == 0:
    raise ValueError("No valid embeddings found matching patterns 'X_dif_nsteps*' or 'X_gcn_nlayers*'")

df_filtered = df.loc[valid_indices].copy()
df_filtered['network_layers'] = network_layers_list
df_filtered['method'] = methods_list

print(f"\n=== Filtered data ===")
print(f"Valid embeddings: {len(df_filtered)}")

# Check for "Metric Type" row and remove it if present
METRIC_TYPE_ROW = "Metric Type"
if METRIC_TYPE_ROW in df_filtered.index:
    print(f"\nRemoving '{METRIC_TYPE_ROW}' row from data")
    df_filtered = df_filtered.drop(METRIC_TYPE_ROW)

# Identify aggregate score columns
aggregate_scores = {}
possible_names = {
    'Total': ['Total', 'total', 'Overall'],
    'Batch correction': ['Batch correction', 'Batch Correction', 'batch_correction', 'BatchCorrection'],
    'Bio conservation': ['Bio conservation', 'Bio Conservation', 'bio_conservation', 'BioConservation']
}

found_scores = {}
for score_name, possible_names_list in possible_names.items():
    for col in df_filtered.columns:
        if col in possible_names_list:
            found_scores[score_name] = col
            break

if len(found_scores) == 0:
    # Try to find columns containing these keywords
    for col in df_filtered.columns:
        col_lower = col.lower()
        if 'total' in col_lower and 'Total' not in found_scores:
            found_scores['Total'] = col
        elif ('batch' in col_lower and 'correction' in col_lower) and 'Batch correction' not in found_scores:
            found_scores['Batch correction'] = col
        elif ('bio' in col_lower and 'conservation' in col_lower) and 'Bio conservation' not in found_scores:
            found_scores['Bio conservation'] = col

if len(found_scores) == 0:
    raise ValueError(f"Could not find aggregate score columns. Available columns: {list(df_filtered.columns)}")

print(f"\n=== Found aggregate scores ===")
for score_name, col_name in found_scores.items():
    print(f"  {score_name}: {col_name}")

# Convert aggregate score columns to numeric
for score_name, col_name in found_scores.items():
    df_filtered[col_name] = pd.to_numeric(df_filtered[col_name], errors='coerce')

# Separate data by method
celldiffusion_df = df_filtered[df_filtered['method'] == 'CellDiffusion'].copy()
gcn_df = df_filtered[df_filtered['method'] == 'GCN'].copy()

# Sort by network_layers
celldiffusion_df = celldiffusion_df.sort_values('network_layers')
gcn_df = gcn_df.sort_values('network_layers')

print(f"\nCellDiffusion data points: {len(celldiffusion_df)}")
print(f"GCN data points: {len(gcn_df)}")

# Create figure with subplots for each aggregate score
n_scores = len(found_scores)
fig, axes = plt.subplots(1, n_scores, figsize=(6 * n_scores, 6))

if n_scores == 1:
    axes = [axes]

# Get all network layers for x-axis
all_layers = sorted(df_filtered['network_layers'].unique())

# Professional color palette
colors = {
    'CellDiffusion': '#4472C4',  # Rich blue
    'GCN': '#C55A11'  # Warm orange-red
}

# Plot each aggregate score
for i, (score_name, col_name) in enumerate(found_scores.items()):
    ax = axes[i]
    
    # Plot CellDiffusion line
    if len(celldiffusion_df) > 0:
        x_values = celldiffusion_df['network_layers'].values
        y_values = pd.to_numeric(celldiffusion_df[col_name], errors='coerce').values
        mask = ~pd.isna(y_values)
        if mask.sum() > 0:
            ax.plot(x_values[mask], y_values[mask], 
                   marker='o', linewidth=2.5, markersize=8, label='CellDiffusion', 
                   color=colors['CellDiffusion'], alpha=0.9, 
                   markeredgecolor='white', markeredgewidth=1.0)
    
    # Plot GCN line
    if len(gcn_df) > 0:
        x_values = gcn_df['network_layers'].values
        y_values = pd.to_numeric(gcn_df[col_name], errors='coerce').values
        mask = ~pd.isna(y_values)
        if mask.sum() > 0:
            ax.plot(x_values[mask], y_values[mask], 
                   marker='s', linewidth=2.5, markersize=8, label='GCN', 
                   color=colors['GCN'], alpha=0.9,
                   markeredgecolor='white', markeredgewidth=1.0)
    
    # Customize subplot
    ax.set_xlabel('Network Layers', fontsize=12, fontweight='bold', labelpad=8)
    ax.set_ylabel('Score', fontsize=12, fontweight='bold', labelpad=8)
    ax.set_title(f'{score_name}', fontsize=14, fontweight='bold', pad=10)
    
    # Improve grid styling
    ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5, axis='both', zorder=0)
    ax.set_axisbelow(True)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)
    
    # Set x-axis ticks
    ax.set_xticks(all_layers)
    ax.set_xticklabels([int(x) for x in all_layers], rotation=0, fontsize=10)
    
    # Set y-axis limits
    if len(df_filtered) > 0:
        values = pd.to_numeric(df_filtered[col_name], errors='coerce').dropna()
        if len(values) > 0:
            y_min = float(values.min())
            y_max = float(values.max())
            y_range = y_max - y_min
            y_padding = y_range * 0.08 if y_range > 0 else 0.1
            ax.set_ylim(bottom=max(0, y_min - y_padding), 
                       top=min(1.0, y_max + y_padding))
        else:
            ax.set_ylim(bottom=0, top=1)
    
    # Add minor ticks
    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))
    
    # Improve tick styling
    ax.tick_params(axis='both', which='major', length=6, width=1.0)
    ax.tick_params(axis='both', which='minor', length=3, width=0.5)
    
    # Add legend to first subplot only
    if i == 0:
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
                           fontsize=10)
        legend.get_frame().set_linewidth(0.8)

# Add overall title
fig.suptitle('SCIB Evaluation (Mesenchymal Lineage): Aggregate Scores vs Network Layers', 
             fontsize=16, fontweight='bold', y=1.02)

plt.tight_layout()

# Save plot
print(f"\n=== Saving plot ===")
print(f"Saving to: {output_pdf}")
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

print("Plot saved successfully!")
print(f"  Compared methods: CellDiffusion, GCN")
print(f"  Aggregate scores: {list(found_scores.keys())}")

print("\n=== Plotting complete! ===")

