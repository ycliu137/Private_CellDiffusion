"""
Plot UMAP visualizations
"""
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib import rcParams

FIGSIZE = (6, 6)
rcParams["figure.figsize"] = FIGSIZE

# Load input data
input_h5ad = snakemake.input.h5ad
output_pdf = snakemake.output.pdf
params = snakemake.params

print(f"Loading data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

print(f"Data shape: {adata.shape}")

# Check which UMAP embeddings are available
print(f"\n=== Available UMAP embeddings ===")
umap_keys_available = []
for key in params.umap_keys:
    if key in adata.obsm:
        umap_keys_available.append(key)
        print(f"  {key}: shape {adata.obsm[key].shape}")
    else:
        print(f"  Warning: {key} not found in adata.obsm")

if len(umap_keys_available) == 0:
    raise ValueError("No UMAP embeddings found in adata.obsm")

# Create figure with subplots
# We'll plot each UMAP with batch, labels, and lineage
# Add an extra row for the shared legend
n_umaps = len(umap_keys_available)
n_cols = 3  # batch, labels, and lineage
n_rows = n_umaps + 1  # +1 for the legend row

fig, axes = plt.subplots(n_rows, n_cols, figsize=(FIGSIZE[0] * n_cols, FIGSIZE[1] * n_rows))
if n_umaps == 1:
    if n_cols == 1:
        axes = [[axes[0]], [axes[1]]]
    else:
        axes = [[axes[0, 0], axes[0, 1], axes[0, 2]], [axes[1, 0], axes[1, 1], axes[1, 2]]]
else:
    if n_cols == 1:
        axes = [[ax] for ax in axes[:-1]] + [[axes[-1]]]
    else:
        axes = axes.reshape(n_rows, n_cols)

# Check if lineage column exists
lineage_key = 'lineage'
has_lineage = lineage_key in adata.obs.columns

# Plot each UMAP
for i, umap_key in enumerate(umap_keys_available):
    # Convert umap_key to friendly display name
    if 'X_umap_dif' in umap_key or 'dif' in umap_key.lower():
        display_name = "CellDiffusion"
    elif 'X_umap_gcn' in umap_key or 'gcn' in umap_key.lower():
        display_name = "GCN"
    else:
        display_name = umap_key  # Fallback to original key if no match
    
    # Temporarily set UMAP for plotting
    adata.obsm['X_umap'] = adata.obsm[umap_key].copy()
    
    # Plot by batch
    sc.pl.umap(
        adata,
        color=params.batch_key,
        ax=axes[i][0],
        show=False,
        frameon=False,
        title=f"{display_name} - Batch",
        legend_loc='none'  # Disable individual legend
    )
    
    # Plot by labels
    if hasattr(params, 'label_key'):
        sc.pl.umap(
            adata,
            color=params.label_key,
            ax=axes[i][1],
            show=False,
            frameon=False,
            title=f"{display_name} - Labels",
            legend_loc='none'  # Disable individual legend
        )
    else:
        axes[i][1].axis('off')
    
    # Plot by lineage
    if has_lineage:
        sc.pl.umap(
            adata,
            color=lineage_key,
            ax=axes[i][2],
            show=False,
            frameon=False,
            title=f"{display_name} - Lineage",
            legend_loc='none'  # Disable individual legend
        )
    else:
        axes[i][2].axis('off')
        print(f"  Warning: '{lineage_key}' column not found in adata.obs, skipping lineage plot")
    
    # Clear temporary UMAP
    del adata.obsm['X_umap']

# Hide all axes in the bottom row (legend row)
axes[n_umaps][0].axis('off')
axes[n_umaps][1].axis('off')
axes[n_umaps][2].axis('off')

# Get legend handles and labels from one of the labels plots
# We'll use the first UMAP's labels plot to get the legend
handles = None
labels = None
if hasattr(params, 'label_key') and len(umap_keys_available) > 0:
    # Temporarily set UMAP for getting legend
    adata.obsm['X_umap'] = adata.obsm[umap_keys_available[0]].copy()
    
    # Create a temporary plot to extract legend handles and labels
    temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
    sc.pl.umap(
        adata,
        color=params.label_key,
        ax=temp_ax,
        show=False,
        frameon=False
    )
    handles, labels = temp_ax.get_legend_handles_labels()
    plt.close(temp_fig)
    
    del adata.obsm['X_umap']
else:
    # If no label_key, try to get legend from batch plot
    if len(umap_keys_available) > 0:
        adata.obsm['X_umap'] = adata.obsm[umap_keys_available[0]].copy()
        temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
        sc.pl.umap(
            adata,
            color=params.batch_key,
            ax=temp_ax,
            show=False,
            frameon=False
        )
        handles, labels = temp_ax.get_legend_handles_labels()
        plt.close(temp_fig)
        
        del adata.obsm['X_umap']

# Apply tight_layout first to adjust subplot positions  
plt.tight_layout()

# Create shared legend spanning all three columns in the bottom row
if handles and labels:
    # Use fig.legend() to place legend in the bottom row, spanning all columns
    # Calculate the y position for the bottom row (approximately at 1/(2*n_rows) from bottom)
    # Using a small offset to place it nicely in the bottom row area
    bottom_y = max(0.01, 0.5 / n_rows)
    fig.legend(handles, labels, loc='lower center', ncol=min(len(labels), 8),
               frameon=False, fontsize=8, bbox_to_anchor=(0.5, bottom_y))

# Save figure
print(f"\n=== Saving UMAP plot ===")
print(f"Saving to: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print("\n=== UMAP plotting complete! ===")

