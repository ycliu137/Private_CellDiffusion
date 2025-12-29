"""
Plot collapse diagnostic: neighbor purity vs network layers
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

print(f"Loading aggregated metrics from: {input_csv}")
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

# Create plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot CellDiffusion line
if len(celldiffusion_df) > 0:
    ax.plot(celldiffusion_df['network_layers'], celldiffusion_df['neighbor_purity'], 
            marker='o', linewidth=2, markersize=8, label='CellDiffusion', 
            color='#1f77b4', alpha=0.8)

# Plot GCN line
if len(gcn_df) > 0:
    ax.plot(gcn_df['network_layers'], gcn_df['neighbor_purity'], 
            marker='s', linewidth=2, markersize=8, label='GCN', 
            color='#ff7f0e', alpha=0.8)

# Customize plot
ax.set_xlabel('Network Layers', fontsize=12, fontweight='bold')
ax.set_ylabel('Neighbor Purity', fontsize=12, fontweight='bold')
ax.set_title('Collapse Diagnostic: Neighbor Purity vs Network Layers', 
             fontsize=14, fontweight='bold')
ax.legend(loc='best', fontsize=11, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--', axis='both')

# Set x-axis ticks to show all network layers values
all_layers = sorted(df['network_layers'].unique())
ax.set_xticks(all_layers)
ax.set_xticklabels([int(x) for x in all_layers])

# Set y-axis limits
if len(df) > 0:
    y_min = df['neighbor_purity'].min()
    y_max = df['neighbor_purity'].max()
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

