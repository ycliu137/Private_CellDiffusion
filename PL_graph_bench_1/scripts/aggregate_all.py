"""
Aggregate X_dif and UMAP embeddings from all graph building method results into a single adata file
Combines aggregate_X_dif.py and aggregate_umap.py
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

print(f"\n=== Aggregating X_dif and UMAP from {len(input_h5ad_files)} files ===")

# Load the first file as the base adata structure
print(f"Loading base adata from: {input_h5ad_files[0]}")
base_adata = sc.read_h5ad(input_h5ad_files[0])
print(f"Base adata shape: {base_adata.shape}")

# Extract X_dif and UMAP embeddings from all files
x_dif_dict = {}
umap_dict = {}

for h5ad_file in input_h5ad_files:
    print(f"\nProcessing: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    filename = Path(h5ad_file).name
    
    # Extract graph method name from filename (format: integrated_{method}.h5ad)
    match = re.search(r'integrated_(.+)\.h5ad', filename)
    if match:
        method_name = match.group(1)
    else:
        print(f"Warning: Could not extract method name from filename {filename}, skipping")
        continue
    
    # Extract X_dif
    if 'X_dif' in adata.obsm:
        key_name = f"X_dif_{method_name}"
        x_dif_dict[key_name] = adata.obsm['X_dif'].copy()
        print(f"  Extracted X_dif shape: {adata.obsm['X_dif'].shape}, stored as '{key_name}'")
    else:
        print(f"  Warning: X_dif not found in {filename}")
    
    # Extract UMAP embeddings
    umap_keys = [k for k in adata.obsm.keys() if 'umap' in k.lower()]
    for umap_key in umap_keys:
        # Create a standardized key name: X_umap_dif_{method}
        new_key = f"X_umap_dif_{method_name}"
        if new_key not in umap_dict:
            umap_dict[new_key] = adata.obsm[umap_key].copy()
            print(f"  Extracted {umap_key} shape: {adata.obsm[umap_key].shape}, stored as '{new_key}'")
        else:
            print(f"  Warning: Key '{new_key}' already exists, skipping")

# Check if we have any embeddings to aggregate
if len(x_dif_dict) == 0:
    raise ValueError(
        "No X_dif embeddings found to aggregate. "
        "Please check that all integration files contain 'X_dif' in adata.obsm."
    )

# Add all X_dif to base adata
print(f"\n=== Adding all X_dif to base adata ===")
for key, x_dif_matrix in x_dif_dict.items():
    base_adata.obsm[key] = x_dif_matrix
    print(f"  Added {key}: shape {x_dif_matrix.shape}")

# Add all UMAP embeddings to base adata
print(f"\n=== Adding all UMAP embeddings to base adata ===")
for key, umap_matrix in umap_dict.items():
    base_adata.obsm[key] = umap_matrix
    print(f"  Added {key}: shape {umap_matrix.shape}")

# Also add metadata columns to track which methods were aggregated
base_adata.uns['aggregated_X_dif'] = {
    'methods': sorted(list(x_dif_dict.keys())),
    'num_representations': len(x_dif_dict)
}

if len(umap_dict) > 0:
    base_adata.uns['aggregated_umap'] = {
        'methods': sorted(list(umap_dict.keys())),
        'num_representations': len(umap_dict)
    }

print(f"\n=== Saving aggregated adata ===")
print(f"  Output file: {output_h5ad}")
print(f"  Total X_dif representations: {len(x_dif_dict)}")
print(f"  Total UMAP representations: {len(umap_dict)}")
print(f"  X_dif methods: {base_adata.uns['aggregated_X_dif']['methods']}")
if len(umap_dict) > 0:
    print(f"  UMAP methods: {base_adata.uns['aggregated_umap']['methods']}")

# Ensure output directory exists
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)

# Save aggregated adata
base_adata.write(output_h5ad)
print(f"  Saved successfully!")

print("\n=== Aggregation complete! ===")

