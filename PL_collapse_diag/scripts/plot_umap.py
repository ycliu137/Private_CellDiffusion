"""
Plot UMAP visualizations for individual integration results
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

# Load input data
input_h5ad = snakemake.input.h5ad
output_pdf = snakemake.output.pdf
params = snakemake.params

print(f"Loading data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

print(f"Data shape: {adata.shape}")
print(f"Method: {params.method}")
print(f"Network layers: {params.network_layers}")

# Determine which UMAP key to use
if params.method == "CellDiffusion":
    umap_key = f'X_umap_dif_nsteps{params.network_layers}'
elif params.method == "GCN":
    umap_key = f'X_umap_gcn_nlayers{params.network_layers}'
else:
    raise ValueError(f"Unknown method: {params.method}")

print(f"Using UMAP key: {umap_key}")

# Check if UMAP exists
if umap_key not in adata.obsm:
    raise ValueError(f"UMAP key '{umap_key}' not found in adata.obsm. Available keys: {[k for k in adata.obsm.keys() if 'umap' in k.lower()]}")

# Create figure with 2 subplots (batch and labels)
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Temporarily set UMAP for plotting
adata.obsm['X_umap'] = adata.obsm[umap_key].copy()

# Plot by batch
sc.pl.umap(
    adata,
    color=params.batch_key,
    ax=axes[0],
    show=False,
    frameon=False,
    title=f"{params.method} (n={params.network_layers}) - Batch"
)

# Plot by labels
sc.pl.umap(
    adata,
    color=params.label_key,
    ax=axes[1],
    show=False,
    frameon=False,
    title=f"{params.method} (n={params.network_layers}) - Labels"
)

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

