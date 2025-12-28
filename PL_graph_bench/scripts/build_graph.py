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
print(f"\n=== Building integration graph using {params.graph_method} ===")

if params.graph_method == 'OMNN-Harmony':
    cd.inte.build_integration_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
elif params.graph_method == 'pure-MNN':
    cd.inte.build_mnn_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
elif params.graph_method == 'Harmony-MNN':
    cd.inte.build_harmony_mnn_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
elif params.graph_method == 'scVI_MNN':
    cd.inte.build_scvi_mnn_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
elif params.graph_method == 'OMNN-scVI':
    cd.inte.build_omnn_scvi_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
else:
    raise ValueError(f"Unknown graph building method: {params.graph_method}")

print(f"Integration graph shape: {adata.uns['integration_edge_index'].shape}")

# Save data with graph
print(f"\n=== Saving data with integration graph ===")
print(f"Saving to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

print(f"\n=== Graph building complete using {params.graph_method}! ===")

