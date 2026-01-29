"""
Plot UMAP visualizations for specific network layers (2, 8, 16) comparing CellDiffusion and GCN
All plots in one figure with shared legends for batch and cell type (labels)
Square subplots for each UMAP visualization - Publication quality
"""
import sys
from pathlib import Path
import re

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Get input file and output from snakemake
try:
    input_h5ad = snakemake.input.h5ad
    output_pdf = snakemake.output.pdf
    params = snakemake.params
except Exception as e:
    print(f"Error getting snakemake inputs/outputs: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print(f"\n=== Loading aggregated UMAP data ===")
print(f"Input file: {input_h5ad}")

# Load aggregated adata
adata = sc.read_h5ad(input_h5ad)
print(f"Data shape: {adata.shape}")

# Find all UMAP keys
print(f"\n=== Finding UMAP embeddings ===")
umap_keys = [k for k in adata.obsm.keys() if 'umap' in k.lower()]
print(f"Found {len(umap_keys)} UMAP embeddings:")
for key in umap_keys:
    print(f"  - {key}")

if len(umap_keys) == 0:
    raise ValueError("No UMAP embeddings found in input data")

# Organize UMAP keys by network layer
layer_data = {}  # {layer_number: {'celldiffusion': umap_key, 'gcn': umap_key}}

for umap_key in umap_keys:
    # Pattern: X_umap_dif_nsteps{N} or X_umap_gcn_nlayers{N}
    if 'X_umap_dif_nsteps' in umap_key:
        match = re.search(r'X_umap_dif_nsteps(\d+)', umap_key)
        if match:
            layer = int(match.group(1))
            if layer not in layer_data:
                layer_data[layer] = {}
            layer_data[layer]['celldiffusion'] = umap_key
            print(f"  Found CellDiffusion layer {layer}: {umap_key}")
    elif 'X_umap_gcn_nlayers' in umap_key:
        match = re.search(r'X_umap_gcn_nlayers(\d+)', umap_key)
        if match:
            layer = int(match.group(1))
            if layer not in layer_data:
                layer_data[layer] = {}
            layer_data[layer]['gcn'] = umap_key
            print(f"  Found GCN layer {layer}: {umap_key}")

# Select specific layers to plot: 2, 8, 16
target_layers = [2, 8, 16]
layers_to_plot = [l for l in target_layers if l in layer_data]

if len(layers_to_plot) == 0:
    print(f"Warning: None of the target layers {target_layers} found in data")
    print(f"Available layers: {sorted(layer_data.keys())}")
    # Use available layers instead
    layers_to_plot = sorted(layer_data.keys())[:3]  # Use first 3 available

print(f"\n=== Selected layers to plot: {layers_to_plot} ===")

# Create figure layout
# 3 rows (one per layer) x 4 columns (CellDiffusion batch, CellDiffusion cell type, GCN batch, GCN cell type)
# Plus 1 row for legends
n_layers = len(layers_to_plot)
n_cols = 4  # CellDiffusion batch, CellDiffusion labels, GCN batch, GCN labels
n_rows = n_layers + 1  # +1 for legend row

# Calculate figure size with square subplots (publication quality)
subplot_size = 6  # Size of each square subplot (in inches)
fig_width = subplot_size * n_cols + 1.5  # Add extra width for spacing
fig_height = subplot_size * n_layers + 2.5  # Height for UMAP rows + legend row

fig = plt.figure(figsize=(fig_width, fig_height))
gs = GridSpec(n_rows, n_cols, figure=fig,
              hspace=0.4, wspace=0.4,  # Spacing between subplots
              left=0.06, right=0.96, top=0.96, bottom=0.08)

# Create axes array for UMAP plots
axes = []
for i in range(n_layers):
    row_axes = []
    for j in range(n_cols):
        ax = fig.add_subplot(gs[i, j])
        ax.set_aspect('equal', adjustable='box')  # Make plots square
        row_axes.append(ax)
    axes.append(row_axes)

# Create axes for legend row
legend_axes = []
for j in range(n_cols):
    ax = fig.add_subplot(gs[n_layers, j])
    legend_axes.append(ax)

# Add main title
fig.suptitle('Network Architecture Comparison: CellDiffusion vs GCN', 
             fontsize=20, fontweight='bold', y=0.98)

# Plot UMAPs for each selected layer
print(f"\n=== Plotting UMAPs ===")
for layer_idx, layer in enumerate(layers_to_plot):
    print(f"\nLayer {layer}:")
    
    # Plot CellDiffusion results (batch and labels)
    if 'celldiffusion' in layer_data[layer]:
        try:
            umap_key_dif = layer_data[layer]['celldiffusion']
            print(f"  Plotting CellDiffusion (batch)...")
            
            if umap_key_dif in adata.obsm:
                # Temporarily set UMAP for plotting
                adata.obsm['X_umap'] = adata.obsm[umap_key_dif].copy()
                
                # Plot CellDiffusion - Batch
                sc.pl.umap(
                    adata,
                    color=params.batch_key,
                    ax=axes[layer_idx][0],
                    show=False,
                    frameon=True,
                    legend_loc='none',
                    title=f"L={layer} | CellDiffusion | Batch" if layer_idx == 0 else f"L={layer}",
                    title_fontsize=14
                )
                axes[layer_idx][0].title.set_fontsize(14)
                axes[layer_idx][0].title.set_fontweight('bold')
                
                # Plot CellDiffusion - Cell Type (Labels)
                print(f"  Plotting CellDiffusion (cell type)...")
                sc.pl.umap(
                    adata,
                    color=params.label_key,
                    ax=axes[layer_idx][1],
                    show=False,
                    frameon=True,
                    legend_loc='none',
                    title=f"CellDiffusion | Cell Type" if layer_idx == 0 else "",
                    title_fontsize=14
                )
                if layer_idx == 0:
                    axes[layer_idx][1].title.set_fontsize(14)
                    axes[layer_idx][1].title.set_fontweight('bold')
                
                # Clean up
                del adata.obsm['X_umap']
            else:
                print(f"    Warning: {umap_key_dif} not found in adata")
                axes[layer_idx][0].text(0.5, 0.5, 'Data not available', ha='center', va='center', fontsize=12)
                axes[layer_idx][1].text(0.5, 0.5, 'Data not available', ha='center', va='center', fontsize=12)
        except Exception as e:
            print(f"    Error plotting CellDiffusion: {e}")
            import traceback
            traceback.print_exc()
            axes[layer_idx][0].text(0.5, 0.5, 'Error', ha='center', va='center', fontsize=12)
            axes[layer_idx][1].text(0.5, 0.5, 'Error', ha='center', va='center', fontsize=12)
    else:
        print(f"  Warning: CellDiffusion data not found for layer {layer}")
        axes[layer_idx][0].text(0.5, 0.5, 'Not available', ha='center', va='center', fontsize=12)
        axes[layer_idx][1].text(0.5, 0.5, 'Not available', ha='center', va='center', fontsize=12)
    
    # Plot GCN results (batch and labels)
    if 'gcn' in layer_data[layer]:
        try:
            umap_key_gcn = layer_data[layer]['gcn']
            print(f"  Plotting GCN (batch)...")
            
            if umap_key_gcn in adata.obsm:
                # Temporarily set UMAP for plotting
                adata.obsm['X_umap'] = adata.obsm[umap_key_gcn].copy()
                
                # Plot GCN - Batch
                sc.pl.umap(
                    adata,
                    color=params.batch_key,
                    ax=axes[layer_idx][2],
                    show=False,
                    frameon=True,
                    legend_loc='none',
                    title=f"GCN | Batch" if layer_idx == 0 else "",
                    title_fontsize=14
                )
                if layer_idx == 0:
                    axes[layer_idx][2].title.set_fontsize(14)
                    axes[layer_idx][2].title.set_fontweight('bold')
                
                # Plot GCN - Cell Type (Labels)
                print(f"  Plotting GCN (cell type)...")
                sc.pl.umap(
                    adata,
                    color=params.label_key,
                    ax=axes[layer_idx][3],
                    show=False,
                    frameon=True,
                    legend_loc='none',
                    title=f"GCN | Cell Type" if layer_idx == 0 else "",
                    title_fontsize=14
                )
                if layer_idx == 0:
                    axes[layer_idx][3].title.set_fontsize(14)
                    axes[layer_idx][3].title.set_fontweight('bold')
                
                # Clean up
                del adata.obsm['X_umap']
            else:
                print(f"    Warning: {umap_key_gcn} not found in adata")
                axes[layer_idx][2].text(0.5, 0.5, 'Data not available', ha='center', va='center', fontsize=12)
                axes[layer_idx][3].text(0.5, 0.5, 'Data not available', ha='center', va='center', fontsize=12)
        except Exception as e:
            print(f"    Error plotting GCN: {e}")
            import traceback
            traceback.print_exc()
            axes[layer_idx][2].text(0.5, 0.5, 'Error', ha='center', va='center', fontsize=12)
            axes[layer_idx][3].text(0.5, 0.5, 'Error', ha='center', va='center', fontsize=12)
    else:
        print(f"  Warning: GCN data not found for layer {layer}")
        axes[layer_idx][2].text(0.5, 0.5, 'Not available', ha='center', va='center', fontsize=12)
        axes[layer_idx][3].text(0.5, 0.5, 'Not available', ha='center', va='center', fontsize=12)

# Extract legends (batch and labels only) from first available layer
print(f"\n=== Creating shared legends ===")
batch_handles, batch_labels = None, None
label_handles, label_labels = None, None

if len(layers_to_plot) > 0:
    # Get first layer with both methods to extract legends
    first_layer = layers_to_plot[0]
    
    # Temporarily set UMAP for legend extraction
    if 'celldiffusion' in layer_data[first_layer]:
        adata.obsm['X_umap'] = adata.obsm[layer_data[first_layer]['celldiffusion']].copy()
    elif 'gcn' in layer_data[first_layer]:
        adata.obsm['X_umap'] = adata.obsm[layer_data[first_layer]['gcn']].copy()
    
    # Extract batch legend
    temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
    sc.pl.umap(
        adata,
        color=params.batch_key,
        ax=temp_ax,
        show=False,
        frameon=False
    )
    batch_handles, batch_labels = temp_ax.get_legend_handles_labels()
    plt.close(temp_fig)
    
    # Extract labels legend
    temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
    sc.pl.umap(
        adata,
        color=params.label_key,
        ax=temp_ax,
        show=False,
        frameon=False
    )
    label_handles, label_labels = temp_ax.get_legend_handles_labels()
    plt.close(temp_fig)
    
    del adata.obsm['X_umap']

# Create shared legends in the legend row
# Columns 0-1: Batch legend
# Columns 2-3: Cell Type legend

if batch_handles and batch_labels:
    legend_axes[0].axis('off')
    legend_axes[1].axis('off')
    # Create batch legend spanning columns 0-1
    batch_legend = fig.legend(batch_handles, batch_labels,
                              loc='lower left',
                              frameon=True, fontsize=12,
                              ncol=5,  # Multiple columns for batch legend
                              bbox_to_anchor=(0.055, 0.01),
                              title='Batch', title_fontsize=13,
                              framealpha=0.98,
                              edgecolor='black', fancybox=True)
else:
    legend_axes[0].axis('off')
    legend_axes[1].axis('off')

if label_handles and label_labels:
    legend_axes[2].axis('off')
    legend_axes[3].axis('off')
    # Create cell type legend spanning columns 2-3
    label_legend = fig.legend(label_handles, label_labels,
                              loc='lower center',
                              frameon=True, fontsize=12,
                              ncol=7,  # Multiple columns for cell type legend
                              bbox_to_anchor=(0.65, 0.01),
                              title='Cell Type', title_fontsize=13,
                              framealpha=0.98,
                              edgecolor='black', fancybox=True)
else:
    legend_axes[2].axis('off')
    legend_axes[3].axis('off')

# Save figure
print(f"\n=== Saving UMAP plot ===")
print(f"Saving to: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print("\n=== UMAP network layers comparison plot complete! ===")
