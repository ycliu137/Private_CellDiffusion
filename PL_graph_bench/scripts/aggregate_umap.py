"""
Aggregate UMAP embeddings from all graph building method results into a single adata file
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

print(f"\n=== Aggregating UMAP embeddings from {len(input_h5ad_files)} files ===")

# Load the first file as the base adata structure
print(f"Loading base adata from: {input_h5ad_files[0]}")
base_adata = sc.read_h5ad(input_h5ad_files[0])
print(f"Base adata shape: {base_adata.shape}")

# Extract UMAP embeddings from all files
umap_dict = {}

for h5ad_file in input_h5ad_files:
    print(f"\nProcessing: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    filename = Path(h5ad_file).name
    
    # Extract graph method name from filename (format: integrated_umap_{method}.h5ad)
    match = re.search(r'integrated_umap_(.+)\.h5ad', filename)
    if match:
        method_name = match.group(1)
    else:
        print(f"Warning: Could not extract method name from filename {filename}, skipping")
        continue
    
    # Look for UMAP keys in the adata
    umap_keys = [k for k in adata.obsm.keys() if 'umap' in k.lower()]
    for umap_key in umap_keys:
        # Create a standardized key name: X_umap_dif_{method}
        new_key = f"X_umap_dif_{method_name}"
        if new_key not in umap_dict:
            umap_dict[new_key] = adata.obsm[umap_key].copy()
            print(f"  Extracted {umap_key} shape: {adata.obsm[umap_key].shape}, stored as '{new_key}'")
        else:
            print(f"  Warning: Key '{new_key}' already exists, skipping")

# Add all UMAP embeddings to base adata
print(f"\n=== Adding all UMAP embeddings to base adata ===")
for key, umap_matrix in umap_dict.items():
    base_adata.obsm[key] = umap_matrix
    print(f"  Added {key}: shape {umap_matrix.shape}")

# Also add a metadata column to track which methods were aggregated
base_adata.uns['aggregated_umap'] = {
    'methods': sorted(list(umap_dict.keys())),
    'num_representations': len(umap_dict)
}

print(f"\n=== Saving aggregated adata ===")
print(f"  Output file: {output_h5ad}")
print(f"  Total UMAP representations: {len(umap_dict)}")
print(f"  Methods: {base_adata.uns['aggregated_umap']['methods']}")

# Ensure output directory exists
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)

# Save aggregated adata
base_adata.write(output_h5ad)
print(f"  Saved successfully!")

print("\n=== Aggregation complete! ===")

