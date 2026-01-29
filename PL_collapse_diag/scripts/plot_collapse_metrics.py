"""
Plot collapse diagnostic metrics: intrinsic_dimension_knn, intrinsic_dimension_at_variance_percentage, variance_explained_by_embedding, neighbor_purity vs network layers
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

# Create output paths for individual plots
output_dir = Path(output_pdf).parent
output_pdf_stem = Path(output_pdf).stem
output_combined = output_pdf
output_intrinsic_knn = output_dir / f"{output_pdf_stem}_intrinsic_dimension_knn.pdf"
output_intrinsic_var = output_dir / f"{output_pdf_stem}_intrinsic_dimension_variance.pdf"
output_variance_exp = output_dir / f"{output_pdf_stem}_variance_explained.pdf"
output_neighbor_pur = output_dir / f"{output_pdf_stem}_neighbor_purity.pdf"

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

# Create figure with four subplots in one row
fig, axes = plt.subplots(1, 4, figsize=(20, 6))

# Define colors
color_celldiffusion = '#1f77b4'
color_gcn = '#ff7f0e'

# Get all network layers for x-axis ticks
all_layers = sorted(df['network_layers'].unique())

# Plot 1: Intrinsic Dimension (KNN)
ax = axes[0]
if len(celldiffusion_df) > 0:
    ax.plot(celldiffusion_df['network_layers'], celldiffusion_df['intrinsic_dimension_knn'], 
            marker='o', linewidth=3, markersize=10, label='CellDiffusion', 
            color=color_celldiffusion, alpha=0.8)
if len(gcn_df) > 0:
    ax.plot(gcn_df['network_layers'], gcn_df['intrinsic_dimension_knn'], 
            marker='s', linewidth=3, markersize=10, label='GCN', 
            color=color_gcn, alpha=0.8)
ax.set_xlabel('Network Layers', fontsize=14, fontweight='bold')
ax.set_ylabel('Intrinsic Dimension (KNN)', fontsize=14, fontweight='bold')
ax.set_title('Intrinsic Dimension (KNN)', fontsize=15, fontweight='bold')
ax.legend(loc='best', fontsize=12, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--', axis='both')
ax.set_xticks(all_layers)
ax.set_xticklabels([int(x) for x in all_layers], fontsize=11)
ax.tick_params(axis='y', labelsize=11)

# Plot 2: Intrinsic Dimension (Variance Percentage)
ax = axes[1]
if len(celldiffusion_df) > 0:
    ax.plot(celldiffusion_df['network_layers'], celldiffusion_df['intrinsic_dimension_variance'], 
            marker='o', linewidth=3, markersize=10, label='CellDiffusion', 
            color=color_celldiffusion, alpha=0.8)
if len(gcn_df) > 0:
    ax.plot(gcn_df['network_layers'], gcn_df['intrinsic_dimension_variance'], 
            marker='s', linewidth=3, markersize=10, label='GCN', 
            color=color_gcn, alpha=0.8)
ax.set_xlabel('Network Layers', fontsize=14, fontweight='bold')
ax.set_ylabel('Intrinsic Dimension (Variance)', fontsize=14, fontweight='bold')
ax.set_title('Intrinsic Dimension (Variance)', fontsize=15, fontweight='bold')
ax.legend(loc='best', fontsize=12, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--', axis='both')
ax.set_xticks(all_layers)
ax.set_xticklabels([int(x) for x in all_layers], fontsize=11)
ax.tick_params(axis='y', labelsize=11)

# Plot 3: Variance Explained
ax = axes[2]
if len(celldiffusion_df) > 0:
    # Filter out NaN values for variance_explained
    celldiffusion_valid = celldiffusion_df[celldiffusion_df['variance_explained'].notna()]
    if len(celldiffusion_valid) > 0:
        ax.plot(celldiffusion_valid['network_layers'], celldiffusion_valid['variance_explained'], 
                marker='o', linewidth=3, markersize=10, label='CellDiffusion', 
                color=color_celldiffusion, alpha=0.8)
if len(gcn_df) > 0:
    # Filter out NaN values for variance_explained
    gcn_valid = gcn_df[gcn_df['variance_explained'].notna()]
    if len(gcn_valid) > 0:
        ax.plot(gcn_valid['network_layers'], gcn_valid['variance_explained'], 
                marker='s', linewidth=3, markersize=10, label='GCN', 
                color=color_gcn, alpha=0.8)
ax.set_xlabel('Network Layers', fontsize=14, fontweight='bold')
ax.set_ylabel('Variance Explained', fontsize=14, fontweight='bold')
ax.set_title('Variance Explained by Embedding', fontsize=15, fontweight='bold')
ax.legend(loc='best', fontsize=12, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--', axis='both')
ax.set_xticks(all_layers)
ax.set_xticklabels([int(x) for x in all_layers], fontsize=11)
ax.tick_params(axis='y', labelsize=11)
# Set y-axis limits for variance explained (should be between 0 and 1)
if len(df[df['variance_explained'].notna()]) > 0:
    y_min = df[df['variance_explained'].notna()]['variance_explained'].min()
    y_max = df[df['variance_explained'].notna()]['variance_explained'].max()
    y_range = y_max - y_min
    ax.set_ylim(bottom=max(0, y_min - 0.05 * y_range), 
                top=min(1.0, y_max + 0.05 * y_range))

# Plot 4: Neighbor Purity
ax = axes[3]
if len(celldiffusion_df) > 0:
    # Filter out NaN values for neighbor_purity
    celldiffusion_valid = celldiffusion_df[celldiffusion_df['neighbor_purity'].notna()]
    if len(celldiffusion_valid) > 0:
        ax.plot(celldiffusion_valid['network_layers'], celldiffusion_valid['neighbor_purity'], 
                marker='o', linewidth=3, markersize=10, label='CellDiffusion', 
                color=color_celldiffusion, alpha=0.8)
if len(gcn_df) > 0:
    # Filter out NaN values for neighbor_purity
    gcn_valid = gcn_df[gcn_df['neighbor_purity'].notna()]
    if len(gcn_valid) > 0:
        ax.plot(gcn_valid['network_layers'], gcn_valid['neighbor_purity'], 
                marker='s', linewidth=3, markersize=10, label='GCN', 
                color=color_gcn, alpha=0.8)
ax.set_xlabel('Network Layers', fontsize=14, fontweight='bold')
ax.set_ylabel('Neighbor Purity', fontsize=14, fontweight='bold')
ax.set_title('Neighbor Purity', fontsize=15, fontweight='bold')
ax.legend(loc='best', fontsize=12, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--', axis='both')
ax.set_xticks(all_layers)
ax.set_xticklabels([int(x) for x in all_layers], fontsize=11)
ax.tick_params(axis='y', labelsize=11)
# Set y-axis limits for neighbor purity (should be between 0 and 1)
if len(df[df['neighbor_purity'].notna()]) > 0:
    y_min = df[df['neighbor_purity'].notna()]['neighbor_purity'].min()
    y_max = df[df['neighbor_purity'].notna()]['neighbor_purity'].max()
    y_range = y_max - y_min
    ax.set_ylim(bottom=max(0, y_min - 0.05 * y_range), 
                top=min(1.0, y_max + 0.05 * y_range))

plt.tight_layout()

# Save combined plot
print(f"\n=== Saving combined plot ===")
print(f"Saving to: {output_combined}")
Path(output_combined).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_combined, dpi=300, bbox_inches='tight')
plt.close()

print("Combined plot saved successfully!")

# Create individual plots for each metric
def create_individual_plot(metric_name, metric_col, ylabel, output_file):
    """Create individual plot for a metric"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot CellDiffusion
    if len(celldiffusion_df) > 0:
        valid_data = celldiffusion_df[celldiffusion_df[metric_col].notna()]
        if len(valid_data) > 0:
            ax.plot(valid_data['network_layers'], valid_data[metric_col], 
                   marker='o', linewidth=4, markersize=12, label='CellDiffusion', 
                   color=color_celldiffusion, alpha=0.8)
    
    # Plot GCN
    if len(gcn_df) > 0:
        valid_data = gcn_df[gcn_df[metric_col].notna()]
        if len(valid_data) > 0:
            ax.plot(valid_data['network_layers'], valid_data[metric_col], 
                   marker='s', linewidth=4, markersize=12, label='GCN', 
                   color=color_gcn, alpha=0.8)
    
    ax.set_xlabel('Network Layers', fontsize=16, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=16, fontweight='bold')
    ax.set_title(metric_name, fontsize=18, fontweight='bold')
    ax.legend(loc='best', fontsize=14, framealpha=0.95, edgecolor='black', fancybox=True)
    ax.grid(True, alpha=0.3, linestyle='--', axis='both')
    ax.set_xticks(all_layers)
    ax.set_xticklabels([int(x) for x in all_layers], fontsize=13)
    ax.tick_params(axis='y', labelsize=13)
    
    # Set y-axis limits for metrics between 0 and 1
    if metric_col in ['variance_explained', 'neighbor_purity']:
        valid_values = df[df[metric_col].notna()][metric_col]
        if len(valid_values) > 0:
            y_min = valid_values.min()
            y_max = valid_values.max()
            y_range = y_max - y_min
            ax.set_ylim(bottom=max(0, y_min - 0.05 * y_range), 
                        top=min(1.0, y_max + 0.05 * y_range))
    
    plt.tight_layout()
    print(f"Saving individual plot to: {output_file}")
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")

# Create individual plots for each metric
create_individual_plot('Intrinsic Dimension (KNN)', 'intrinsic_dimension_knn', 
                       'Intrinsic Dimension (KNN)', output_intrinsic_knn)
create_individual_plot('Intrinsic Dimension (Variance)', 'intrinsic_dimension_variance', 
                       'Intrinsic Dimension (Variance)', output_intrinsic_var)
create_individual_plot('Variance Explained by Embedding', 'variance_explained', 
                       'Variance Explained', output_variance_exp)
create_individual_plot('Neighbor Purity', 'neighbor_purity', 
                       'Neighbor Purity', output_neighbor_pur)

print("\n=== All plots saved successfully ===")

