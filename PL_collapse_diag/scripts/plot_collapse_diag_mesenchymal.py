"""
Plot collapse diagnostic for Mesenchymal lineage: neighbor purity vs network layers
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set publication-quality style parameters
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans', 'Liberation Sans']
rcParams['font.size'] = 11
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 14
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 10
rcParams['axes.linewidth'] = 1.0
rcParams['grid.linewidth'] = 0.5
rcParams['lines.linewidth'] = 2.0
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

# Get input/output from snakemake
input_csv = snakemake.input.csv
output_pdf = snakemake.output.pdf

print(f"\n=== Loading aggregated metrics (Mesenchymal lineage) ===")
print(f"Input file: {input_csv}")
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

# Professional color palette
colors = {
    'CellDiffusion': '#4472C4',  # Rich blue
    'GCN': '#C55A11'  # Warm orange-red
}

# Plot CellDiffusion line
if len(celldiffusion_df) > 0:
    ax.plot(celldiffusion_df['network_layers'], celldiffusion_df['neighbor_purity'], 
            marker='o', linewidth=2.5, markersize=8, label='CellDiffusion', 
            color=colors['CellDiffusion'], alpha=0.9,
            markeredgecolor='white', markeredgewidth=1.0)

# Plot GCN line
if len(gcn_df) > 0:
    ax.plot(gcn_df['network_layers'], gcn_df['neighbor_purity'], 
            marker='s', linewidth=2.5, markersize=8, label='GCN', 
            color=colors['GCN'], alpha=0.9,
            markeredgecolor='white', markeredgewidth=1.0)

# Customize plot
ax.set_xlabel('Network Layers', fontsize=12, fontweight='bold', labelpad=8)
ax.set_ylabel('Neighbor Purity', fontsize=12, fontweight='bold', labelpad=8)
ax.set_title('Collapse Diagnostic (Mesenchymal Lineage): Neighbor Purity vs Network Layers', 
             fontsize=14, fontweight='bold', pad=10)

# Improve grid styling
ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5, axis='both', zorder=0)
ax.set_axisbelow(True)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.0)
ax.spines['bottom'].set_linewidth(1.0)

# Legend
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

# Set x-axis ticks to show all network layers values
all_layers = sorted(df['network_layers'].unique())
ax.set_xticks(all_layers)
ax.set_xticklabels([int(x) for x in all_layers], fontsize=10)

# Add minor ticks
ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))

# Improve tick styling
ax.tick_params(axis='both', which='major', length=6, width=1.0)
ax.tick_params(axis='both', which='minor', length=3, width=0.5)

# Set y-axis limits
if len(df) > 0:
    y_min = df['neighbor_purity'].min()
    y_max = df['neighbor_purity'].max()
    y_range = y_max - y_min
    y_padding = y_range * 0.08 if y_range > 0 else 0.05
    ax.set_ylim(bottom=max(0, y_min - y_padding), 
                top=min(1.0, y_max + y_padding))

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

