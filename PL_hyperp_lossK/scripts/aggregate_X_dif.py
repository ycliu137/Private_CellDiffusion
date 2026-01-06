"""
Aggregate X_dif from all loss_k integration outputs into a single adata file
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

print(f"\n=== Aggregating X_dif from {len(input_h5ad_files)} files ===")

# Load the first file as the base adata structure
print(f"Loading base adata from: {input_h5ad_files[0]}")
base_adata = sc.read_h5ad(input_h5ad_files[0])
print(f"Base adata shape: {base_adata.shape}")

# Extract loss_k values from filenames and store X_dif
x_dif_dict = {}

for h5ad_file in input_h5ad_files:
    print(f"\nProcessing: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    
    # Extract loss_k value from filename (format: integrated_lossk{loss_k}.h5ad)
    filename = Path(h5ad_file).name
    
    # Try to extract loss_k from filename
    match = re.search(r'lossk(\d+)', filename)
    if match:
        loss_k = int(match.group(1))
    else:
        print(f"Warning: Could not extract loss_k from filename {filename}, skipping")
        continue
    
    # Check if X_dif exists
    if 'X_dif' not in adata.obsm:
        print(f"Warning: X_dif not found in {filename}, skipping")
        continue
    
    # Store X_dif with key indicating loss_k value
    key_name = f"X_dif_lossk{loss_k}"
    x_dif_dict[key_name] = adata.obsm['X_dif'].copy()
    print(f"  Extracted X_dif shape: {adata.obsm['X_dif'].shape}, stored as '{key_name}'")

# Check if we have any X_dif embeddings to aggregate
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

# Also add a metadata column to track which loss_k values were aggregated
base_adata.uns['aggregated_X_dif'] = {
    'loss_k_values': sorted([int(k.split('lossk')[1]) for k in x_dif_dict.keys()]),
    'num_representations': len(x_dif_dict)
}

print(f"\n=== Saving aggregated adata ===")
print(f"  Output file: {output_h5ad}")
print(f"  Total X_dif representations: {len(x_dif_dict)}")
print(f"  loss_k values: {base_adata.uns['aggregated_X_dif']['loss_k_values']}")

# Ensure output directory exists
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)

# Save aggregated adata
base_adata.write(output_h5ad)
print(f"  Saved successfully!")

print("\n=== Aggregation complete! ===")
