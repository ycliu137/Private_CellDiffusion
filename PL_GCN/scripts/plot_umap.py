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
# We'll plot each UMAP with both batch and labels
n_umaps = len(umap_keys_available)
n_cols = 2  # batch and labels
n_rows = n_umaps

fig, axes = plt.subplots(n_rows, n_cols, figsize=(FIGSIZE[0] * n_cols, FIGSIZE[1] * n_rows))
if n_umaps == 1:
    if n_cols == 1:
        axes = [[axes]]
    else:
        axes = [[axes[0], axes[1]]]
else:
    if n_cols == 1:
        axes = [[ax] for ax in axes]
    else:
        axes = axes.reshape(n_rows, n_cols)

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
        title=f"{display_name} - Batch"
    )
    
    # Plot by labels
    if hasattr(params, 'label_key'):
        sc.pl.umap(
            adata,
            color=params.label_key,
            ax=axes[i][1],
            show=False,
            frameon=False,
            title=f"{display_name} - Labels"
        )
    else:
        axes[i][1].axis('off')
    
    # Clear temporary UMAP
    del adata.obsm['X_umap']

plt.tight_layout()

# Save figure
print(f"\n=== Saving UMAP plot ===")
print(f"Saving to: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print("\n=== UMAP plotting complete! ===")

