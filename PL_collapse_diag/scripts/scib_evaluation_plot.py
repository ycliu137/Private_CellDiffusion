"""
Plot SCIB evaluation: three aggregate scores vs network layers for CellDiffusion and GCN
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import re

# Get input/output from snakemake
input_csv = snakemake.input.csv
output_pdf = snakemake.output.pdf

print(f"Loading SCIB results table from: {input_csv}")
df = pd.read_csv(input_csv, index_col=0)

print(f"\n=== Data summary ===")
print(f"Results table shape: {df.shape}")
print(f"\nFirst few rows:")
print(df.head())
print(f"\nColumn names: {list(df.columns)}")
print(f"\nIndex (Embedding names): {list(df.index)}")

# Extract network layers and method from embedding names
# Format: X_dif_nsteps{N} or X_gcn_nlayers{N}
network_layers_list = []
methods_list = []
valid_indices = []

for idx in df.index:
    # Try to match CellDiffusion pattern: X_dif_nsteps{N}
    match_dif = re.search(r'X_dif_nsteps(\d+)', idx)
    # Try to match GCN pattern: X_gcn_nlayers{N}
    match_gcn = re.search(r'X_gcn_nlayers(\d+)', idx)
    
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
df_filtered = df.loc[valid_indices].copy()
df_filtered['network_layers'] = network_layers_list
df_filtered['method'] = methods_list

print(f"\n=== Filtered data ===")
print(f"Valid embeddings: {len(df_filtered)}")
print(df_filtered[['network_layers', 'method']].head())

# Identify aggregate score columns
# Look for: Total, Batch correction, Bio conservation
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
        if 'total' in col_lower and 'total' not in found_scores:
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

# Convert aggregate score columns to numeric (handle any string values)
print(f"\n=== Converting score columns to numeric ===")
for score_name, col_name in found_scores.items():
    print(f"  Converting {col_name} to numeric...")
    df_filtered[col_name] = pd.to_numeric(df_filtered[col_name], errors='coerce')
    # Check for any NaN values after conversion
    nan_count = df_filtered[col_name].isna().sum()
    if nan_count > 0:
        print(f"    Warning: {nan_count} NaN values found in {col_name} after conversion")

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

# Plot each aggregate score
for i, (score_name, col_name) in enumerate(found_scores.items()):
    ax = axes[i]
    
    # Plot CellDiffusion line
    if len(celldiffusion_df) > 0:
        # Ensure values are numeric
        x_values = celldiffusion_df['network_layers'].values
        y_values = pd.to_numeric(celldiffusion_df[col_name], errors='coerce').values
        # Filter out NaN values
        mask = ~pd.isna(y_values)
        if mask.sum() > 0:
            ax.plot(x_values[mask], y_values[mask], 
                    marker='o', linewidth=2.5, markersize=10, label='CellDiffusion', 
                    color='#1f77b4', alpha=0.8, markeredgecolor='black', markeredgewidth=1.2)
    
    # Plot GCN line
    if len(gcn_df) > 0:
        # Ensure values are numeric
        x_values = gcn_df['network_layers'].values
        y_values = pd.to_numeric(gcn_df[col_name], errors='coerce').values
        # Filter out NaN values
        mask = ~pd.isna(y_values)
        if mask.sum() > 0:
            ax.plot(x_values[mask], y_values[mask], 
                    marker='s', linewidth=2.5, markersize=10, label='GCN', 
                    color='#ff7f0e', alpha=0.8, markeredgecolor='black', markeredgewidth=1.2)
    
    # Customize subplot
    ax.set_xlabel('Network Layers', fontsize=12, fontweight='bold')
    ax.set_ylabel('Score', fontsize=12, fontweight='bold')
    ax.set_title(f'{score_name}', fontsize=13, fontweight='bold')
    ax.legend(loc='best', fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle='--', axis='both')
    
    # Set x-axis ticks
    ax.set_xticks(all_layers)
    ax.set_xticklabels([int(x) for x in all_layers], rotation=0)
    
    # Set y-axis limits
    if len(df_filtered) > 0:
        # Ensure column is numeric and drop NaN values
        values = pd.to_numeric(df_filtered[col_name], errors='coerce').dropna()
        if len(values) > 0:
            y_min = float(values.min())
            y_max = float(values.max())
            y_range = y_max - y_min
            ax.set_ylim(bottom=max(0, y_min - 0.05 * y_range), 
                        top=min(1.0, y_max + 0.05 * y_range))
        else:
            print(f"    Warning: No valid numeric values found for {col_name}, using default y-axis limits")
            ax.set_ylim(bottom=0, top=1)

# Add overall title
fig.suptitle('SCIB Evaluation: Aggregate Scores vs Network Layers', 
             fontsize=16, fontweight='bold', y=1.02)

plt.tight_layout()

# Save plot
print(f"\n=== Saving plot ===")
print(f"Saving to: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print("Plot saved successfully!")

