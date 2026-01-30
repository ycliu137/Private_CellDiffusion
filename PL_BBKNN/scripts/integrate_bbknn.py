"""
BBKNN integration
"""
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import bbknn
import pandas as pd
import matplotlib.pyplot as plt

input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
output_umap = snakemake.output.umap
params = snakemake.params

print(f"\n=== Loading data ===")
print(f"  Input: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)
print(f"  Data shape: {adata.shape}")

# Preprocessing
print(f"\n=== Step 0: Preprocessing ===")
print(f"  Normalized data: {params.normalized_data}")

if not params.normalized_data:
    print(f"  Filtering genes with min_cells={params.min_cells}")
    sc.pp.filter_genes(adata, min_cells=params.min_cells)
    
    print(f"  Normalizing with target_sum={params.target_sum}")
    sc.pp.normalize_total(adata, target_sum=params.target_sum)
    sc.pp.log1p(adata)
    print(f"  Normalization complete")
else:
    print("  Data already normalized, skipping preprocessing")

# HVG selection
print(f"  Finding highly variable genes: n_top_genes={params.n_top_genes}")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=params.n_top_genes,
    min_mean=params.min_mean,
    max_mean=params.max_mean,
    min_disp=params.min_disp
)

adata.raw = adata.copy()
adata = adata[:, adata.var.highly_variable]
print(f"  After HVG selection: {adata.shape}")

# Scaling
print(f"\n=== Step 1: Scaling ===")
sc.pp.scale(adata, max_value=10)
print(f"  Scaling complete")

# PCA for neighbor computation
print(f"\n=== Step 2: Computing PCA ===")
sc.tl.pca(adata, svd_solver='arpack', n_comps=params.n_pcs)
print(f"  PCA complete")

# BBKNN integration
print(f"\n=== Step 3: BBKNN integration ===")
print(f"  Batch key: {params.batch_key}")
print(f"  Computation: {params.computation}")

bbknn.bbknn(
    adata,
    batch_key=params.batch_key,
    computation=params.computation
)
print(f"  BBKNN integration complete")

# UMAP
print(f"\n=== Step 4: Computing UMAP ===")
sc.tl.umap(adata, min_dist=params.umap_min_dist)
adata.obsm['X_umap_bbknn'] = adata.obsm['X_umap'].copy()
print(f"  UMAP complete! Saved to adata.obsm['X_umap_bbknn']")

# Leiden clustering
print(f"\n=== Step 5: Leiden clustering ===")
print(f"  Resolution: {params.leiden_resolution}")
sc.tl.leiden(adata, resolution=params.leiden_resolution, key_added='leiden_bbknn')
print(f"  Leiden clustering complete! Saved to adata.obs['leiden_bbknn']")

# Save integrated data
print(f"\n=== Saving integrated data ===")
print(f"  Output: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

# Plot UMAP
print(f"\n=== Creating UMAP plot ===")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

sc.pl.umap(adata, color=params.batch_key, ax=axes[0], show=False)
axes[0].set_title('BBKNN - Batch')

sc.pl.umap(adata, color=params.label_key, ax=axes[1], show=False)
axes[1].set_title('BBKNN - Labels')

plt.tight_layout()
print(f"  Saving UMAP plot to: {output_umap}")
Path(output_umap).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_umap, dpi=150, bbox_inches='tight')
plt.close()

print("\n=== BBKNN integration complete! ===")
