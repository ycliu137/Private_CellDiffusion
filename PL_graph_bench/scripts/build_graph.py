"""
Build integration graph using different graph building methods
"""
import sys
from pathlib import Path
import pandas as pd

# Add project root to path first, before importing celldiffusion
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import celldiffusion as cd
import torch
import time
import random

# Load input data
input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
params = snakemake.params

print(f"Loading encoded data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

print(f"Data shape: {adata.shape}")
print(f"Using graph building method: {params.graph_method}")
print(f"Graph method type: {type(params.graph_method)}")
print(f"Graph method repr: {repr(params.graph_method)}")

# Normalize method name to handle any potential parsing issues
method_name = str(params.graph_method).strip()

# Check if adata.X exists (required by SCVI-based methods)
print(f"\n=== Checking data availability ===")
print(f"adata.X exists: {adata.X is not None}")
if adata.X is not None:
    print(f"adata.X shape: {adata.X.shape}")
    print(f"adata.X type: {type(adata.X)}")
    if hasattr(adata.X, 'toarray'):
        print(f"adata.X is sparse matrix")
    print(f"adata.X dtype: {adata.X.dtype if hasattr(adata.X, 'dtype') else 'N/A'}")
else:
    print("WARNING: adata.X is None!")

# Check if adata.raw exists (backup for SCVI if needed)
if hasattr(adata, 'raw') and adata.raw is not None:
    print(f"adata.raw exists: True")
    print(f"adata.raw.X shape: {adata.raw.X.shape}")
else:
    print(f"adata.raw exists: False")

# CUDA synchronization: Wait for GPU to be available
if params.device == "cuda":
    print("\n=== Waiting for GPU availability ===")
    if torch.cuda.is_available():
        initial_delay = random.uniform(0, 2)
        time.sleep(initial_delay)
        wait_interval = 5
        waited = 0
        while True:
            try:
                test_tensor = torch.zeros(1, device='cuda')
                torch.cuda.synchronize()
                del test_tensor
                torch.cuda.empty_cache()
                print(f"GPU is available (waited {waited:.1f}s)")
                break
            except RuntimeError as e:
                if "busy" in str(e).lower() or "unavailable" in str(e).lower():
                    time.sleep(wait_interval)
                    waited += wait_interval
                else:
                    raise

# Create node_batch_mt (required by some graph building functions)
print("\n=== Creating node_batch_mt ===")
if 'node_batch_mt' not in adata.obsm:
    adata.obsm['node_batch_mt'] = pd.get_dummies(adata.obs[params.batch_key]).to_numpy()
    print("Created node_batch_mt")
else:
    print("node_batch_mt already exists")

# Build integration loss adjacency (required by build_integration_graph)
print("\n=== Building integration loss adjacency ===")
cd.inte.build_integration_loss_adj(
    adata,
    use_rep='X_fae',
    k=params.k,
    device=params.device
)

# Build graph using the specified method
print(f"\n=== Building integration graph using {method_name} ===")

if method_name == 'OMNN-Harmony':
    cd.inte.build_integration_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
elif method_name == 'pure-MNN':
    cd.inte.build_mnn_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
elif method_name == 'Harmony-MNN':
    cd.inte.build_harmony_mnn_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
elif method_name == 'scVI-MNN' or method_name == 'scVI_MNN':
    # SCVI needs adata.X to exist (it will use it for training)
    # Note: SCVI prefers raw count data, but can work with normalized/log-transformed data
    if adata.X is None:
        raise ValueError("adata.X is None! SCVI requires adata.X to train the model. "
                        "Please ensure the preprocessing step preserves adata.X.")
    print(f"adata.X is available for SCVI training (shape: {adata.X.shape})")
    print("Note: SCVI may warn if data is not raw counts, but it should still work.")
    cd.inte.build_scvi_mnn_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
elif method_name == 'OMNN-scVI':
    # SCVI needs adata.X to exist (it will use it for training)
    if adata.X is None:
        raise ValueError("adata.X is None! SCVI requires adata.X to train the model. "
                        "Please ensure the preprocessing step preserves adata.X.")
    print(f"adata.X is available for SCVI training (shape: {adata.X.shape})")
    print("Note: SCVI may warn if data is not raw counts, but it should still work.")
    cd.inte.build_omnn_scvi_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
else:
    raise ValueError(f"Unknown graph building method: {method_name} (original: {params.graph_method}, type: {type(params.graph_method)})")

print(f"Integration graph shape: {adata.uns['integration_edge_index'].shape}")

# Save data with graph
print(f"\n=== Saving data with integration graph ===")
print(f"Saving to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

print(f"\n=== Graph building complete using {method_name}! ===")

