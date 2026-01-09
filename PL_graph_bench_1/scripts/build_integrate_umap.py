"""
Build graph, integrate, and compute UMAP in one step to reduce file I/O
Combines build_graph.py, integrate.py, and compute_umap.py
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

# Normalize method name
method_name = str(params.graph_method).strip()

# Check if adata.X exists (required by SCVI-based methods)
print(f"\n=== Checking data availability ===")
print(f"adata.X exists: {adata.X is not None}")
if adata.X is not None:
    print(f"adata.X shape: {adata.X.shape}")
else:
    print("WARNING: adata.X is None!")

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

# ===== STEP 1: Build Graph =====
print("\n=== STEP 1: Building integration graph ===")

# Create node_batch_mt (required by some graph building functions)
if 'node_batch_mt' not in adata.obsm:
    adata.obsm['node_batch_mt'] = pd.get_dummies(adata.obs[params.batch_key]).to_numpy()
    print("Created node_batch_mt")
else:
    print("node_batch_mt already exists")

# Build integration loss adjacency (required by build_integration_graph)
cd.inte.build_integration_loss_adj(
    adata,
    use_rep='X_fae',
    k=params.k,
    device=params.device
)

# Build graph using the specified method
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
    if adata.X is None:
        raise ValueError("adata.X is None! SCVI requires adata.X to train the model.")
    print(f"adata.X is available for SCVI training (shape: {adata.X.shape})")
    cd.inte.build_scvi_mnn_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
elif method_name == 'OMNN-scVI':
    if adata.X is None:
        raise ValueError("adata.X is None! SCVI requires adata.X to train the model.")
    print(f"adata.X is available for SCVI training (shape: {adata.X.shape})")
    cd.inte.build_omnn_scvi_graph(
        adata,
        batch_key=params.batch_key,
        use_rep='X_fae',
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device
    )
else:
    raise ValueError(f"Unknown graph building method: {method_name}")

print(f"Integration graph shape: {adata.uns['integration_edge_index'].shape}")

# ===== STEP 2: Run CellDiffusion Integration =====
print("\n=== STEP 2: Running CellDiffusion integration ===")
print(f"  Max epochs: {params.max_epoch}")
print(f"  Learning rate: {params.lr}")
print(f"  Num features: {params.num_features_diffusion}")
print(f"  Num heads: {params.num_heads_diffusion}")
print(f"  Num steps: {params.num_steps_diffusion}")
print(f"  Time increment: {params.time_increment_diffusion}")
print(f"  Device: {params.device}")

cd.inte.integration_diffusion(
    adata,
    use_rep='X_fae',
    save_key='X_dif',  # Save as X_dif
    max_epoch=params.max_epoch,
    lr=params.lr,
    num_features_diffusion=params.num_features_diffusion,
    num_heads_diffusion=params.num_heads_diffusion,
    num_steps_diffusion=params.num_steps_diffusion,
    time_increment_diffusion=params.time_increment_diffusion,
    device=params.device
)

print(f"CellDiffusion integration complete! Embeddings saved in adata.obsm['X_dif']")

# ===== STEP 3: Compute UMAP =====
print("\n=== STEP 3: Computing UMAP ===")
print(f"  Use rep: {params.use_rep}")
print(f"  UMAP key: {params.umap_key}")
print(f"  n_neighbors: {params.n_neighbors}")
print(f"  n_pcs: {params.n_pcs}")

if params.use_rep not in adata.obsm:
    print(f"  Warning: {params.use_rep} not found in adata.obsm, skipping UMAP")
else:
    # Compute neighbors and UMAP
    sc.pp.neighbors(
        adata,
        use_rep=params.use_rep,
        n_neighbors=params.n_neighbors,
        n_pcs=params.n_pcs
    )
    sc.tl.umap(adata)
    
    # Save UMAP to specified key
    if 'X_umap' in adata.obsm:
        adata.obsm[params.umap_key] = adata.obsm['X_umap'].copy()
        print(f"  Saved UMAP to adata.obsm['{params.umap_key}']")
        
        # Clean up temporary neighbor graph and umap
        if 'neighbors' in adata.uns:
            del adata.uns['neighbors']
        if 'umap' in adata.uns:
            del adata.uns['umap']
        if 'X_umap' in adata.obsm:
            del adata.obsm['X_umap']
    else:
        print(f"  Warning: UMAP not computed for {params.use_rep}")

# Save final data with graph, X_dif, and UMAP
print(f"\n=== Saving final data ===")
print(f"Saving to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

print(f"\n=== Complete pipeline for {method_name} finished! ===")
print(f"  - Graph built and stored in adata.uns['integration_edge_index']")
print(f"  - X_dif stored in adata.obsm['X_dif']")
if params.umap_key in adata.obsm:
    print(f"  - UMAP stored in adata.obsm['{params.umap_key}']")

