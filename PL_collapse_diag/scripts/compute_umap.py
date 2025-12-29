"""
Compute UMAP for individual integration results
"""
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc

# Load input data
input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
params = snakemake.params

print(f"Loading data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

print(f"Data shape: {adata.shape}")
print(f"Method: {params.method}")
print(f"Network layers: {params.network_layers}")

# Determine which representation to use
if params.method == "CellDiffusion":
    use_rep = 'X_dif'
    umap_key = f'X_umap_dif_nsteps{params.network_layers}'
elif params.method == "GCN":
    use_rep = 'X_gcn'
    umap_key = f'X_umap_gcn_nlayers{params.network_layers}'
else:
    raise ValueError(f"Unknown method: {params.method}")

print(f"Using representation: {use_rep}")
print(f"UMAP will be saved as: {umap_key}")

# Check if representation exists
if use_rep not in adata.obsm:
    raise ValueError(f"Representation '{use_rep}' not found in adata.obsm. Available keys: {list(adata.obsm.keys())}")

# Compute neighbors and UMAP
print(f"\n=== Computing UMAP ===")
print(f"  Use rep: {use_rep}")
print(f"  n_neighbors: {params.n_neighbors}")
print(f"  n_pcs: {params.n_pcs}")

sc.pp.neighbors(
    adata,
    use_rep=use_rep,
    n_neighbors=params.n_neighbors,
    n_pcs=params.n_pcs
)
sc.tl.umap(adata)

# Save UMAP to specified key
if 'X_umap' in adata.obsm:
    adata.obsm[umap_key] = adata.obsm['X_umap'].copy()
    print(f"  Saved UMAP to adata.obsm['{umap_key}']")
    
    # Clean up temporary neighbor graph and umap
    if 'neighbors' in adata.uns:
        del adata.uns['neighbors']
    if 'umap' in adata.uns:
        del adata.uns['umap']
    if 'X_umap' in adata.obsm:
        del adata.obsm['X_umap']
else:
    raise ValueError(f"UMAP not computed for {use_rep}")

# Save data with UMAP embedding
print(f"\n=== Saving data with UMAP embedding ===")
print(f"Saving to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

print("\n=== UMAP computation complete! ===")

