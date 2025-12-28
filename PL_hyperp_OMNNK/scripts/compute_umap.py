"""
Compute UMAP for embeddings
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

# Check which embeddings are available
print(f"\n=== Available embeddings ===")
for key in adata.obsm.keys():
    if 'X_' in key:
        print(f"  {key}: shape {adata.obsm[key].shape}")

# Compute UMAP for each embedding
print(f"\n=== Computing UMAP ===")
for use_rep, umap_key in zip(params.use_reps, params.umap_keys):
    if use_rep not in adata.obsm:
        print(f"  Warning: {use_rep} not found in adata.obsm, skipping")
        continue
    
    print(f"\n  Computing UMAP for {use_rep} -> {umap_key}")
    print(f"    n_neighbors: {params.n_neighbors}")
    print(f"    n_pcs: {params.n_pcs}")
    
    # Compute neighbors and UMAP
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
        print(f"    Saved UMAP to adata.obsm['{umap_key}']")
        
        # Clean up temporary neighbor graph and umap
        if 'neighbors' in adata.uns:
            del adata.uns['neighbors']
        if 'umap' in adata.uns:
            del adata.uns['umap']
        if 'X_umap' in adata.obsm:
            del adata.obsm['X_umap']
    else:
        print(f"    Warning: UMAP not computed for {use_rep}")

# Save data with UMAP embeddings
print(f"\n=== Saving data with UMAP embeddings ===")
print(f"Saving to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

print("\n=== UMAP computation complete! ===")

