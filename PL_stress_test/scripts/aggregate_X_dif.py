"""
Aggregate X_dif from all integration outputs into a single adata file
"""
import sys
from pathlib import Path
import scanpy as sc

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Get input files and output file from snakemake
input_h5ad_files = snakemake.input.h5ad_files
output_h5ad = snakemake.output.h5ad

print(f"\n=== Aggregating X_dif from {len(input_h5ad_files)} files ===")

# Load the first file as the base adata structure
print(f"Loading base adata from: {input_h5ad_files[0]}")
base_adata = sc.read_h5ad(input_h5ad_files[0])
print(f"Base adata shape: {base_adata.shape}")

# Extract k_add values from filenames and store X_dif
x_dif_dict = {}

for h5ad_file in input_h5ad_files:
    print(f"\nProcessing: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    
    # Extract k_add value from filename (format: scDiffusion_integration_kadd{k_add}.h5ad)
    filename = Path(h5ad_file).name
    if 'kadd' in filename:
        # Extract k_add value from filename
        try:
            k_add_str = filename.split('kadd')[1].split('.h5ad')[0]
            k_add = int(k_add_str)
        except (IndexError, ValueError):
            print(f"Warning: Could not extract k_add from filename {filename}, skipping")
            continue
    else:
        print(f"Warning: Could not find kadd in filename {filename}, skipping")
        continue
    
    # Check if X_dif exists
    if 'X_dif' not in adata.obsm:
        print(f"Warning: X_dif not found in {filename}, skipping")
        continue
    
    # Store X_dif with key indicating k_add value
    key_name = f"X_dif_kadd{k_add}"
    x_dif_dict[key_name] = adata.obsm['X_dif'].copy()
    print(f"  Extracted X_dif shape: {adata.obsm['X_dif'].shape}, stored as '{key_name}'")

# Add all X_dif to base adata
print(f"\n=== Adding all X_dif to base adata ===")
for key, x_dif_matrix in x_dif_dict.items():
    base_adata.obsm[key] = x_dif_matrix
    print(f"  Added {key}: shape {x_dif_matrix.shape}")

# Also add a metadata column to track which k_add values were aggregated
base_adata.uns['aggregated_X_dif'] = {
    'k_add_values': sorted([int(k.split('kadd')[1]) for k in x_dif_dict.keys()]),
    'num_representations': len(x_dif_dict)
}

print(f"\n=== Saving aggregated adata ===")
print(f"  Output file: {output_h5ad}")
print(f"  Total X_dif representations: {len(x_dif_dict)}")
print(f"  k_add values: {base_adata.uns['aggregated_X_dif']['k_add_values']}")

# Ensure output directory exists
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)

# Save aggregated adata
base_adata.write(output_h5ad)
print(f"  Saved successfully!")

print("\n=== Aggregation complete! ===")

