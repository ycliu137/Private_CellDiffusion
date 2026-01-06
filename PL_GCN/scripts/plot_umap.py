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

# Extract legend handles and labels for batch, labels, and lineage
# We'll use the first UMAP to extract legends
batch_handles, batch_labels = None, None
label_handles, label_labels = None, None
lineage_handles, lineage_labels = None, None

if len(umap_keys_available) > 0:
    # Temporarily set UMAP for getting legends
    adata.obsm['X_umap'] = adata.obsm[umap_keys_available[0]].copy()
    
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
    if hasattr(params, 'label_key'):
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
    
    # Extract lineage legend
    if has_lineage:
        temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
        sc.pl.umap(
            adata,
            color=lineage_key,
            ax=temp_ax,
            show=False,
            frameon=False
        )
        lineage_handles, lineage_labels = temp_ax.get_legend_handles_labels()
        plt.close(temp_fig)
    
    del adata.obsm['X_umap']

# Create three legends in the bottom row, each in its own column
# Column 0: Batch legend
if batch_handles and batch_labels:
    axes[n_umaps][0].axis('off')
    batch_legend = axes[n_umaps][0].legend(batch_handles, batch_labels, loc='center',
                                           frameon=False, fontsize=8, 
                                           ncol=min(len(batch_labels), 6))
    # Add title for batch legend
    axes[n_umaps][0].text(0.5, 0.95, 'Batch', transform=axes[n_umaps][0].transAxes,
                          ha='center', va='top', fontsize=10, fontweight='bold')
else:
    axes[n_umaps][0].axis('off')

# Column 1: Labels legend
if label_handles and label_labels:
    axes[n_umaps][1].axis('off')
    label_legend = axes[n_umaps][1].legend(label_handles, label_labels, loc='center',
                                           frameon=False, fontsize=8,
                                           ncol=min(len(label_labels), 6))
    # Add title for labels legend
    axes[n_umaps][1].text(0.5, 0.95, 'Labels', transform=axes[n_umaps][1].transAxes,
                          ha='center', va='top', fontsize=10, fontweight='bold')
else:
    axes[n_umaps][1].axis('off')

# Column 2: Lineage legend
if lineage_handles and lineage_labels:
    axes[n_umaps][2].axis('off')
    lineage_legend = axes[n_umaps][2].legend(lineage_handles, lineage_labels, loc='center',
                                             frameon=False, fontsize=8,
                                             ncol=min(len(lineage_labels), 6))
    # Add title for lineage legend
    axes[n_umaps][2].text(0.5, 0.95, 'Lineage', transform=axes[n_umaps][2].transAxes,
                          ha='center', va='top', fontsize=10, fontweight='bold')
else:
    axes[n_umaps][2].axis('off')

# Apply tight_layout to adjust subplot positions after adding legends
plt.tight_layout()

# Save figure
print(f"\n=== Saving UMAP plot ===")
print(f"Saving to: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print("\n=== UMAP plotting complete! ===")

