"""
Preprocess single-cell data
"""
import sys
from pathlib import Path
import scanpy as sc
import anndata as ad

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Load input data
input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
params = snakemake.params

print(f"Loading data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

print(f"Original data shape: {adata.shape}")

# Filter genes
print(f"Filtering genes with min_cells={params.min_cells}")
sc.pp.filter_genes(adata, min_cells=params.min_cells)

# Normalize and log transform
print(f"Normalizing with target_sum={params.target_sum}")
sc.pp.normalize_total(adata, target_sum=params.target_sum)
sc.pp.log1p(adata)

# Find highly variable genes
print(f"Finding highly variable genes: n_top_genes={params.n_top_genes}")
sc.pp.highly_variable_genes(
    adata, 
    n_top_genes=params.n_top_genes,
    min_mean=params.min_mean,
    max_mean=params.max_mean,
    min_disp=params.min_disp
)

# Store raw data and subset to highly variable genes
adata.raw = adata.copy()
adata = adata[:, adata.var.highly_variable]

print(f"Final data shape: {adata.shape}")

# Save preprocessed data
print(f"Saving preprocessed data to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

print("Preprocessing complete!")

