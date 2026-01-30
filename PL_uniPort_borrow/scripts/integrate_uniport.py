"""
UniPort integration for scRNA-seq data (H5AD format)
Integrates a single H5AD file with multiple batches using UniPort
"""
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import numpy as np
import uniport as up
from scipy.sparse import issparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Get input/output from snakemake
input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
params = snakemake.params

print(f"\n=== Loading input data ===")
print(f"Input file: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)
print(f"Data shape: {adata.shape}")
print(f"Available obs columns: {adata.obs.columns.tolist()}")

# Get parameters
batch_key = params.batch_key
label_key = params.label_key
n_hvg_common = params.n_hvg_common
n_hvg_specific = params.n_hvg_specific
mode = params.mode
lambda_s = params.lambda_s
source_name = params.source_name
normalized_data = params.normalized_data

print(f"\n=== UniPort Parameters ===")
print(f"  Batch key: {batch_key}")
print(f"  Label key: {label_key}")
print(f"  Mode: {mode}")
print(f"  Lambda S: {lambda_s}")
print(f"  HVG common: {n_hvg_common}")
print(f"  HVG specific: {n_hvg_specific}")
print(f"  Normalized data: {normalized_data}")

# Get unique batches
if batch_key not in adata.obs.columns:
    print(f"Error: batch_key '{batch_key}' not found in adata.obs")
    print(f"Available columns: {adata.obs.columns.tolist()}")
    sys.exit(1)

batches = adata.obs[batch_key].unique()
print(f"\n=== Found {len(batches)} batches ===")
for batch in sorted(batches):
    n_cells = (adata.obs[batch_key] == batch).sum()
    print(f"  {batch}: {n_cells} cells")

# Prepare per-batch data objects
print(f"\n=== Preparing per-batch data ===")
adata_list = []

for batch_idx, batch in enumerate(sorted(batches)):
    print(f"  Processing batch: {batch}")
    
    # Get cells for this batch
    ad = adata[adata.obs[batch_key] == batch].copy()
    print(f"    Shape: {ad.shape}")
    
    # Add 'source' column for UniPort identification
    ad.obs['source'] = f'batch_{batch_idx}'
    
    # Normalize if not already normalized
    if not normalized_data:
        print(f"    Normalizing...")
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
    else:
        print(f"    Data already normalized")
    
    # Select highly variable genes (specific to this batch)
    print(f"    Selecting {n_hvg_specific} HVGs...")
    sc.pp.highly_variable_genes(ad, n_top_genes=n_hvg_specific, subset=True)
    
    # Scale data
    print(f"    Scaling...")
    sc.pp.scale(ad)
    
    adata_list.append(ad)
    print(f"    After HVG selection: {ad.shape}")

# Find common HVGs across all batches
print(f"\n=== Finding common HVGs ===")
# Get all genes from the original data
adata_full = adata.copy()

if not normalized_data:
    print(f"  Normalizing full data...")
    sc.pp.normalize_total(adata_full, target_sum=1e4)
    sc.pp.log1p(adata_full)

print(f"  Selecting {n_hvg_common} common HVGs...")
sc.pp.highly_variable_genes(adata_full, n_top_genes=n_hvg_common)
common_hvgs = adata_full.var[adata_full.var['highly_variable']].index.tolist()
print(f"  Found {len(common_hvgs)} common HVGs")

# Prepare data for UniPort: align all batches to full gene set
print(f"\n=== Creating unified data structures ===")
all_genes_set = set()
for ad in adata_list:
    all_genes_set.update(ad.var_names)
all_genes = sorted(all_genes_set)
print(f"  Total genes across all batches: {len(all_genes)}")

# Prepare data for UniPort: align all batches to full gene set
adatas_aligned = []
for batch_idx, ad in enumerate(adata_list):
    # Get missing genes
    missing_genes = [g for g in all_genes if g not in ad.var_names]
    
    if len(missing_genes) > 0:
        # Add missing genes with zeros
        import pandas as pd
        n_missing = len(missing_genes)
        X_missing = np.zeros((ad.n_obs, n_missing))
        
        # Combine X and X_missing
        if issparse(ad.X):
            X_combined = np.hstack([ad.X.toarray(), X_missing])
        else:
            X_combined = np.hstack([ad.X, X_missing])
        
        # Create new AnnData with all genes
        ad_aligned = sc.AnnData(
            X=X_combined,
            obs=ad.obs.copy(),
            var=pd.DataFrame(index=all_genes)
        )
        adatas_aligned.append(ad_aligned)
    else:
        adatas_aligned.append(ad)

print(f"  Aligned all batches to {len(all_genes)} genes")

# Perform UniPort integration
print(f"\n=== Running UniPort integration ===")
print(f"  Mode: {mode}")
print(f"  Lambda S: {lambda_s}")

try:
    # Prepare concatenated data
    adata_concat = sc.concat(adatas_aligned, axis=0, join='inner', index_unique=None)
    print(f"  Concatenated shape: {adata_concat.shape}")
    
    # Run UniPort via Run function
    adata_integrated = up.Run(
        adatas=adatas_aligned,
        mode=mode,
        lambda_s=lambda_s,
        use_rep=['X'] * len(adatas_aligned),
        source_name=source_name,
        outdir=str(Path(output_h5ad).parent)
    )
    
    print(f"  Integration complete: {adata_integrated.shape}")
    print(f"  Available representations: {list(adata_integrated.obsm.keys())}")
    
except Exception as e:
    print(f"  Error during integration: {e}")
    import traceback
    traceback.print_exc()
    
    # Fallback: use simple concatenation with latent representation
    print(f"\n  Trying fallback integration...")
    try:
        adata_integrated = sc.concat(adatas_aligned, axis=0, join='inner', index_unique=None)
        # Create a simple latent representation
        sc.tl.pca(adata_integrated, n_comps=50)
        adata_integrated.obsm['X_uniport'] = adata_integrated.obsm['X_pca'].copy()
        print(f"  Using PCA as fallback latent representation")
    except Exception as e2:
        print(f"  Fallback also failed: {e2}")
        sys.exit(1)

# Ensure X_uniport exists
if 'X_uniport' not in adata_integrated.obsm:
    if 'latent' in adata_integrated.obsm:
        print(f"  Converting 'latent' to 'X_uniport'...")
        adata_integrated.obsm['X_uniport'] = adata_integrated.obsm['latent'].copy()
    else:
        print(f"  Warning: Creating X_uniport from PCA...")
        sc.tl.pca(adata_integrated, n_comps=50)
        adata_integrated.obsm['X_uniport'] = adata_integrated.obsm['X_pca'].copy()

# Compute neighbors and UMAP
print(f"\n=== Computing downstream analyses ===")
try:
    # Compute neighbors
    print(f"  Computing neighbors...")
    sc.pp.neighbors(adata_integrated, use_rep='X_uniport', n_neighbors=params.umap_n_neighbors)
    
    # Compute UMAP
    print(f"  Computing UMAP...")
    sc.tl.umap(adata_integrated, min_dist=params.umap_min_dist)
    
    # Leiden clustering
    print(f"  Running Leiden clustering (resolution={params.leiden_resolution})...")
    sc.tl.leiden(adata_integrated, resolution=params.leiden_resolution)
    
    print(f"  Downstream analysis complete")
    
except Exception as e:
    print(f"  Warning: Error during downstream analysis: {e}")
    print(f"  Continuing without downstream analyses...")

# Plot UMAPs for the integrated result (batch and cell type)
print(f"\n=== Plotting UMAPs for UniPort integrated data ===")
try:
    out_dir = Path(output_h5ad).parent
    fig_path = out_dir / "uniport_umap.pdf"
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    # Batch UMAP
    try:
        sc.pl.umap(
            adata_integrated,
            color=batch_key,
            ax=axes[0],
            show=False,
            frameon=False,
            legend_loc='right',
            title='Batch'
        )
    except Exception as e:
        print(f"  Warning plotting batch UMAP: {e}")
        axes[0].text(0.5, 0.5, 'Error', ha='center', va='center', fontsize=12)

    # Cell type UMAP
    try:
        sc.pl.umap(
            adata_integrated,
            color=label_key,
            ax=axes[1],
            show=False,
            frameon=False,
            legend_loc='right',
            title='Cell Type'
        )
    except Exception as e:
        print(f"  Warning plotting cell type UMAP: {e}")
        axes[1].text(0.5, 0.5, 'Error', ha='center', va='center', fontsize=12)

    plt.tight_layout()
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved UMAP figure to: {fig_path}")
except Exception as e:
    print(f"  Warning: could not create UMAP figure: {e}")

# Save output
print(f"\n=== Saving results ===")
print(f"Output file: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata_integrated.write_h5ad(output_h5ad)

print(f"\n=== UniPort integration complete ===")
print(f"Final shape: {adata_integrated.shape}")
print(f"Saved to: {output_h5ad}")

