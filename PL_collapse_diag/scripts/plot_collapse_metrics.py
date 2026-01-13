"""
Plot collapse diagnostic metrics: intrinsic_dimension_knn, intrinsic_dimension_at_variance_percentage, variance_explained_by_embedding vs network layers
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# Get input/output from snakemake
input_csv = snakemake.input.csv
output_pdf = snakemake.output.pdf

print(f"Loading aggregated collapse metrics from: {input_csv}")
df = pd.read_csv(input_csv)

print(f"\n=== Data summary ===")
print(df.head(10))

# Separate data by method
celldiffusion_df = df[df['method'] == 'CellDiffusion'].copy()
gcn_df = df[df['method'] == 'GCN'].copy()

print(f"\nCellDiffusion data points: {len(celldiffusion_df)}")
print(f"GCN data points: {len(gcn_df)}")

# Sort by network_layers
celldiffusion_df = celldiffusion_df.sort_values('network_layers')
gcn_df = gcn_df.sort_values('network_layers')

# Create figure with three subplots
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Define colors
color_celldiffusion = '#1f77b4'
color_gcn = '#ff7f0e'

# Get all network layers for x-axis ticks
all_layers = sorted(df['network_layers'].unique())

# Plot 1: Intrinsic Dimension (KNN)
ax = axes[0]
if len(celldiffusion_df) > 0:
    ax.plot(celldiffusion_df['network_layers'], celldiffusion_df['intrinsic_dimension_knn'], 
            marker='o', linewidth=2, markersize=8, label='CellDiffusion', 
            color=color_celldiffusion, alpha=0.8)
if len(gcn_df) > 0:
    ax.plot(gcn_df['network_layers'], gcn_df['intrinsic_dimension_knn'], 
            marker='s', linewidth=2, markersize=8, label='GCN', 
            color=color_gcn, alpha=0.8)
ax.set_xlabel('Network Layers', fontsize=12, fontweight='bold')
ax.set_ylabel('Intrinsic Dimension (KNN)', fontsize=12, fontweight='bold')
ax.set_title('Intrinsic Dimension (KNN)', fontsize=13, fontweight='bold')
ax.legend(loc='best', fontsize=10, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--', axis='both')
ax.set_xticks(all_layers)
ax.set_xticklabels([int(x) for x in all_layers])

# Plot 2: Intrinsic Dimension (Variance Percentage)
ax = axes[1]
if len(celldiffusion_df) > 0:
    ax.plot(celldiffusion_df['network_layers'], celldiffusion_df['intrinsic_dimension_variance'], 
            marker='o', linewidth=2, markersize=8, label='CellDiffusion', 
            color=color_celldiffusion, alpha=0.8)
if len(gcn_df) > 0:
    ax.plot(gcn_df['network_layers'], gcn_df['intrinsic_dimension_variance'], 
            marker='s', linewidth=2, markersize=8, label='GCN', 
            color=color_gcn, alpha=0.8)
ax.set_xlabel('Network Layers', fontsize=12, fontweight='bold')
ax.set_ylabel('Intrinsic Dimension (Variance)', fontsize=12, fontweight='bold')
ax.set_title('Intrinsic Dimension (Variance)', fontsize=13, fontweight='bold')
ax.legend(loc='best', fontsize=10, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--', axis='both')
ax.set_xticks(all_layers)
ax.set_xticklabels([int(x) for x in all_layers])

# Plot 3: Variance Explained
ax = axes[2]
if len(celldiffusion_df) > 0:
    # Filter out NaN values for variance_explained
    celldiffusion_valid = celldiffusion_df[celldiffusion_df['variance_explained'].notna()]
    if len(celldiffusion_valid) > 0:
        ax.plot(celldiffusion_valid['network_layers'], celldiffusion_valid['variance_explained'], 
                marker='o', linewidth=2, markersize=8, label='CellDiffusion', 
                color=color_celldiffusion, alpha=0.8)
if len(gcn_df) > 0:
    # Filter out NaN values for variance_explained
    gcn_valid = gcn_df[gcn_df['variance_explained'].notna()]
    if len(gcn_valid) > 0:
        ax.plot(gcn_valid['network_layers'], gcn_valid['variance_explained'], 
                marker='s', linewidth=2, markersize=8, label='GCN', 
                color=color_gcn, alpha=0.8)
ax.set_xlabel('Network Layers', fontsize=12, fontweight='bold')
ax.set_ylabel('Variance Explained', fontsize=12, fontweight='bold')
ax.set_title('Variance Explained by Embedding', fontsize=13, fontweight='bold')
ax.legend(loc='best', fontsize=10, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--', axis='both')
ax.set_xticks(all_layers)
ax.set_xticklabels([int(x) for x in all_layers])
# Set y-axis limits for variance explained (should be between 0 and 1)
if len(df[df['variance_explained'].notna()]) > 0:
    y_min = df[df['variance_explained'].notna()]['variance_explained'].min()
    y_max = df[df['variance_explained'].notna()]['variance_explained'].max()
    y_range = y_max - y_min
    ax.set_ylim(bottom=max(0, y_min - 0.05 * y_range), 
                top=min(1.0, y_max + 0.05 * y_range))

plt.tight_layout()

# Save plot
print(f"\n=== Saving plot ===")
print(f"Saving to: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print("Plot saved successfully!")

