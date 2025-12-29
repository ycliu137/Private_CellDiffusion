"""
Aggregate all X_dif and X_gcn embeddings from different network layers into a single adata file
"""
import sys
from pathlib import Path
import scanpy as sc
import re

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Get input files and output file from snakemake
input_h5ad_files = snakemake.input.h5ad_files
output_h5ad = snakemake.output.h5ad

print(f"\n=== Aggregating embeddings from {len(input_h5ad_files)} files ===")

# Load the first file as the base adata structure
print(f"Loading base adata from: {input_h5ad_files[0]}")
base_adata = sc.read_h5ad(input_h5ad_files[0])
print(f"Base adata shape: {base_adata.shape}")

# Extract embeddings from all files
for h5ad_file in input_h5ad_files:
    print(f"\nProcessing: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    
    # Extract method and network layers from filename
    filename = Path(h5ad_file).name
    
    # Pattern: celldiffusion_nsteps{N}.h5ad or gcn_nlayers{N}.h5ad
    if 'celldiffusion_nsteps' in filename:
        match = re.search(r'celldiffusion_nsteps(\d+)\.h5ad', filename)
        if match:
            method = 'CellDiffusion'
            network_layers = match.group(1)
            embedding_key = 'X_dif'
            new_key = f'X_dif_nsteps{network_layers}'
        else:
            print(f"Warning: Could not extract network layers from filename {filename}, skipping")
            continue
    elif 'gcn_nlayers' in filename:
        match = re.search(r'gcn_nlayers(\d+)\.h5ad', filename)
        if match:
            method = 'GCN'
            network_layers = match.group(1)
            embedding_key = 'X_gcn'
            new_key = f'X_gcn_nlayers{network_layers}'
        else:
            print(f"Warning: Could not extract network layers from filename {filename}, skipping")
            continue
    else:
        print(f"Warning: Unknown file pattern: {filename}, skipping")
        continue
    
    # Check if embedding exists
    if embedding_key not in adata.obsm:
        print(f"Warning: {embedding_key} not found in {filename}, skipping")
        continue
    
    # Store embedding with new key
    base_adata.obsm[new_key] = adata.obsm[embedding_key].copy()
    print(f"  Extracted {embedding_key} shape: {adata.obsm[embedding_key].shape}, stored as '{new_key}'")

# Also aggregate UMAP embeddings if available
print(f"\n=== Aggregating UMAP embeddings ===")
for h5ad_file in input_h5ad_files:
    print(f"Checking UMAP in: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    filename = Path(h5ad_file).name
    
    # Look for UMAP keys
    umap_keys = [k for k in adata.obsm.keys() if 'umap' in k.lower()]
    for umap_key in umap_keys:
        if umap_key not in base_adata.obsm:
            base_adata.obsm[umap_key] = adata.obsm[umap_key].copy()
            print(f"  Added UMAP: {umap_key}")

# Save aggregated adata
print(f"\n=== Saving aggregated adata ===")
print(f"  Output file: {output_h5ad}")
print(f"  Total embedding keys: {len([k for k in base_adata.obsm.keys() if 'X_dif' in k or 'X_gcn' in k])}")
print(f"  Total UMAP keys: {len([k for k in base_adata.obsm.keys() if 'umap' in k.lower()])}")

# Ensure output directory exists
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)

# Save aggregated adata
base_adata.write(output_h5ad)
print(f"  Saved successfully!")

print("\n=== Aggregation complete! ===")

