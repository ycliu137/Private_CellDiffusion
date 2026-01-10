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
input_h5ad_files = snakemake.input.h5ad_files if 'h5ad_files' in snakemake.input else []
input_gcn_files = snakemake.input.gcn_files if 'gcn_files' in snakemake.input else []
all_input_files = list(input_h5ad_files) + list(input_gcn_files)

print(f"\n=== Aggregating X_dif, X_gcn and UMAP from {len(all_input_files)} files ===")
print(f"  CellDiffusion files: {len(input_h5ad_files)}")
print(f"  GCN files: {len(input_gcn_files)}")

# Load the first file as the base adata structure
if len(all_input_files) == 0:
    raise ValueError("No input files provided!")
    
print(f"Loading base adata from: {all_input_files[0]}")
base_adata = sc.read_h5ad(all_input_files[0])
print(f"Base adata shape: {base_adata.shape}")

# Extract X_dif, X_gcn and UMAP embeddings from all files
x_dif_dict = {}
x_gcn_dict = {}
umap_dict = {}

for h5ad_file in all_input_files:
    print(f"\nProcessing: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    filename = Path(h5ad_file).name
    
    # Extract graph method name from filename (format: integrated_{method}.h5ad or gcn_{method}_nlayers{N}.h5ad)
    method_match = re.search(r'integrated_(.+)\.h5ad', filename)
    gcn_match = re.search(r'gcn_(.+)_nlayers(\d+)\.h5ad', filename)
    
    if method_match:
        # CellDiffusion file
        method_name = method_match.group(1)
        
        # Extract X_dif
        if 'X_dif' in adata.obsm:
            key_name = f"X_dif_{method_name}"
            x_dif_dict[key_name] = adata.obsm['X_dif'].copy()
            print(f"  Extracted X_dif shape: {adata.obsm['X_dif'].shape}, stored as '{key_name}'")
        else:
            print(f"  Warning: X_dif not found in {filename}")
        
        # Extract UMAP embeddings for CellDiffusion
        umap_keys = [k for k in adata.obsm.keys() if 'umap' in k.lower() and 'gcn' not in k.lower()]
        for umap_key in umap_keys:
            # Create a standardized key name: X_umap_dif_{method}
            new_key = f"X_umap_dif_{method_name}"
            if new_key not in umap_dict:
                umap_dict[new_key] = adata.obsm[umap_key].copy()
                print(f"  Extracted {umap_key} shape: {adata.obsm[umap_key].shape}, stored as '{new_key}'")
            else:
                print(f"  Warning: Key '{new_key}' already exists, skipping")
                
    elif gcn_match:
        # GCN file
        method_name = gcn_match.group(1)
        num_layers = gcn_match.group(2)
        
        # Extract X_gcn
        if 'X_gcn' in adata.obsm:
            key_name = f"X_gcn_{method_name}_nlayers{num_layers}"
            x_gcn_dict[key_name] = adata.obsm['X_gcn'].copy()
            print(f"  Extracted X_gcn shape: {adata.obsm['X_gcn'].shape}, stored as '{key_name}'")
        else:
            print(f"  Warning: X_gcn not found in {filename}")
        
        # Extract UMAP embeddings for GCN
        umap_keys = [k for k in adata.obsm.keys() if 'umap' in k.lower() and 'gcn' in k.lower()]
        for umap_key in umap_keys:
            # Use the key directly if it follows the pattern
            if umap_key not in umap_dict:
                umap_dict[umap_key] = adata.obsm[umap_key].copy()
                print(f"  Extracted {umap_key} shape: {adata.obsm[umap_key].shape}, stored as '{umap_key}'")
            else:
                print(f"  Warning: Key '{umap_key}' already exists, skipping")
    else:
        print(f"Warning: Could not extract method name from filename {filename}, skipping")
        continue

# Check if we have any embeddings to aggregate
if len(x_dif_dict) == 0 and len(x_gcn_dict) == 0:
    raise ValueError(
        "No embeddings found to aggregate. "
        "Please check that all integration files contain 'X_dif' or 'X_gcn' in adata.obsm."
    )

# Add all X_dif to base adata
if len(x_dif_dict) > 0:
    print(f"\n=== Adding all X_dif to base adata ===")
    for key, x_dif_matrix in x_dif_dict.items():
        base_adata.obsm[key] = x_dif_matrix
        print(f"  Added {key}: shape {x_dif_matrix.shape}")

# Add all X_gcn to base adata
if len(x_gcn_dict) > 0:
    print(f"\n=== Adding all X_gcn to base adata ===")
    for key, x_gcn_matrix in x_gcn_dict.items():
        base_adata.obsm[key] = x_gcn_matrix
        print(f"  Added {key}: shape {x_gcn_matrix.shape}")

# Add all UMAP embeddings to base adata
print(f"\n=== Adding all UMAP embeddings to base adata ===")
for key, umap_matrix in umap_dict.items():
    base_adata.obsm[key] = umap_matrix
    print(f"  Added {key}: shape {umap_matrix.shape}")

# Also add metadata columns to track which methods were aggregated
if len(x_dif_dict) > 0:
    base_adata.uns['aggregated_X_dif'] = {
        'methods': sorted(list(x_dif_dict.keys())),
        'num_representations': len(x_dif_dict)
    }

if len(x_gcn_dict) > 0:
    base_adata.uns['aggregated_X_gcn'] = {
        'methods': sorted(list(x_gcn_dict.keys())),
        'num_representations': len(x_gcn_dict)
    }

if len(umap_dict) > 0:
    base_adata.uns['aggregated_umap'] = {
        'methods': sorted(list(umap_dict.keys())),
        'num_representations': len(umap_dict)
    }

print(f"\n=== Saving aggregated adata ===")
print(f"  Output file: {output_h5ad}")
print(f"  Total X_dif representations: {len(x_dif_dict)}")
print(f"  Total X_gcn representations: {len(x_gcn_dict)}")
print(f"  Total UMAP representations: {len(umap_dict)}")
if len(x_dif_dict) > 0:
    print(f"  X_dif methods: {base_adata.uns['aggregated_X_dif']['methods']}")
if len(x_gcn_dict) > 0:
    print(f"  X_gcn methods: {base_adata.uns['aggregated_X_gcn']['methods']}")
if len(umap_dict) > 0:
    print(f"  UMAP methods: {base_adata.uns['aggregated_umap']['methods']}")

# Ensure output directory exists
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)

# Save aggregated adata
base_adata.write(output_h5ad)
print(f"  Saved successfully!")

print("\n=== Aggregation complete! ===")
