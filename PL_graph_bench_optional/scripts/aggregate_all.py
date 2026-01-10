"""
Aggregate X_dif and UMAP embeddings from all graph building method results into a single adata file
Combines aggregate_X_dif.py and aggregate_umap.py
"""
import sys
from pathlib import Path
import scanpy as sc
import re
import traceback

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Wrap entire script in try-except for better error reporting
try:
    # Get input files and output file from snakemake
    # Snakemake named inputs are accessed as attributes, not dictionary keys
    input_h5ad_files = getattr(snakemake.input, 'h5ad_files', [])
    input_gcn_files = getattr(snakemake.input, 'gcn_files', [])
    output_h5ad = snakemake.output.h5ad

    # Ensure inputs are lists
    if not isinstance(input_h5ad_files, list):
        input_h5ad_files = [input_h5ad_files] if input_h5ad_files else []
    if not isinstance(input_gcn_files, list):
        input_gcn_files = [input_gcn_files] if input_gcn_files else []

    all_input_files = list(input_h5ad_files) + list(input_gcn_files)

    print(f"\n=== Aggregating X_dif, X_gcn and UMAP from {len(all_input_files)} files ===")
    print(f"  CellDiffusion files: {len(input_h5ad_files)}")
    print(f"  GCN files: {len(input_gcn_files)}")
    print(f"  All input files: {all_input_files}")

    # Load the first file as the base adata structure
    if len(all_input_files) == 0:
        raise ValueError("No input files provided!")

    # Prefer using first CellDiffusion file as base, otherwise use first file
    base_file = None
    for f in all_input_files:
        filename = Path(f).name
        if re.search(r'integrated_(.+)\.h5ad', filename):
            base_file = f
            break

    if base_file is None:
        base_file = all_input_files[0]
        
    print(f"Loading base adata from: {base_file}")
    try:
        base_adata = sc.read_h5ad(base_file)
        print(f"Base adata shape: {base_adata.shape}")
    except Exception as e:
        print(f"ERROR: Failed to load base file {base_file}")
        print(f"Error: {e}")
        raise

    # Extract X_dif, X_gcn and UMAP embeddings from all files
    x_dif_dict = {}
    x_gcn_dict = {}
    umap_dict = {}

    for h5ad_file in all_input_files:
        print(f"\nProcessing: {h5ad_file}")
        filename = Path(h5ad_file).name
        
        # Skip if this is the base file and we already loaded it
        if h5ad_file == base_file:
            print(f"  This is the base file, using already loaded adata")
            adata = base_adata
        else:
            try:
                adata = sc.read_h5ad(h5ad_file)
            except Exception as e:
                print(f"ERROR: Failed to load {h5ad_file}")
                print(f"Error: {e}")
                raise
        
        print(f"  File adata shape: {adata.shape}")
        print(f"  Available obsm keys: {list(adata.obsm.keys())}")
        
        # Extract graph method name from filename (format: integrated_{method}.h5ad or gcn_{method}.h5ad)
        method_match = re.search(r'integrated_(.+)\.h5ad', filename)
        gcn_match = re.search(r'gcn_(.+)\.h5ad', filename)
        
        if method_match:
            print(f"  Method match (CellDiffusion): {method_match.group(1)}")
        if gcn_match:
            print(f"  Method match (GCN): {gcn_match.group(1)}")
        
        if method_match:
            # CellDiffusion file
            method_name = method_match.group(1)
            
            # Extract X_dif
            if 'X_dif' in adata.obsm:
                # Label as CellDiffusion
                key_name = f"X_dif_{method_name}"
                x_dif_dict[key_name] = adata.obsm['X_dif'].copy()
                print(f"  Extracted X_dif shape: {adata.obsm['X_dif'].shape}, stored as '{key_name}' (CellDiffusion)")
            else:
                print(f"  Warning: X_dif not found in {filename}")
            
            # Extract UMAP embeddings for CellDiffusion
            umap_keys = [k for k in adata.obsm.keys() if 'umap' in k.lower() and 'gcn' not in k.lower()]
            for umap_key in umap_keys:
                # Create a standardized key name: X_umap_dif_{method}
                new_key = f"X_umap_dif_{method_name}"
                if new_key not in umap_dict:
                    umap_dict[new_key] = adata.obsm[umap_key].copy()
                    print(f"  Extracted {umap_key} shape: {adata.obsm[umap_key].shape}, stored as '{new_key}' (CellDiffusion)")
                else:
                    print(f"  Warning: Key '{new_key}' already exists, skipping")
                    
        elif gcn_match:
            # GCN file
            method_name = gcn_match.group(1)
            
            # Extract X_gcn
            if 'X_gcn' in adata.obsm:
                # Label as GCN
                key_name = f"X_gcn_{method_name}"
                x_gcn_dict[key_name] = adata.obsm['X_gcn'].copy()
                print(f"  Extracted X_gcn shape: {adata.obsm['X_gcn'].shape}, stored as '{key_name}' (GCN)")
            else:
                print(f"  Warning: X_gcn not found in {filename}")
            
            # Extract UMAP embeddings for GCN
            umap_keys = [k for k in adata.obsm.keys() if 'umap' in k.lower() and 'gcn' in k.lower()]
            for umap_key in umap_keys:
                # Standardize key name: X_umap_gcn_{method}
                new_key = f"X_umap_gcn_{method_name}"
                if new_key not in umap_dict:
                    umap_dict[new_key] = adata.obsm[umap_key].copy()
                    print(f"  Extracted {umap_key} shape: {adata.obsm[umap_key].shape}, stored as '{new_key}' (GCN)")
                else:
                    print(f"  Warning: Key '{new_key}' already exists, skipping")
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
    try:
        base_adata.write(output_h5ad)
        print(f"  Saved successfully!")
    except Exception as e:
        print(f"ERROR: Failed to save aggregated file to {output_h5ad}")
        print(f"Error: {e}")
        raise

    print("\n=== Aggregation complete! ===")

except Exception as e:
    print(f"\n=== FATAL ERROR ===")
    print(f"Error type: {type(e).__name__}")
    print(f"Error message: {str(e)}")
    print(f"\n=== Full traceback ===")
    traceback.print_exc()
    sys.exit(1)
