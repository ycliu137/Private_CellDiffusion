"""
Run complete CellDiffusion integration pipeline
"""
import sys
from pathlib import Path

# Add project root to path first
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import celldiffusion as cd
import torch
import time
import random
import pandas as pd

# Load input data
input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
params = snakemake.params

print(f"\n=== Loading data ===")
print(f"  Input: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)
print(f"  Data shape: {adata.shape}")

# Step 0: Preprocessing (if not normalized)
print(f"\n=== Step 0: Preprocessing ===")
print(f"  Normalized data: {params.normalized_data}")

if not params.normalized_data:
    # Filter genes
    print(f"  Filtering genes with min_cells={params.min_cells}")
    sc.pp.filter_genes(adata, min_cells=params.min_cells)
    print(f"  After filtering: {adata.shape}")
    
    # Normalize and log transform
    print(f"  Normalizing with target_sum={params.target_sum}")
    sc.pp.normalize_total(adata, target_sum=params.target_sum)
    sc.pp.log1p(adata)
    print(f"  Normalization complete")
else:
    print("  Data already normalized, skipping preprocessing")

# Find highly variable genes
print(f"  Finding highly variable genes: n_top_genes={params.n_top_genes}")
sc.pp.highly_variable_genes(
    adata, 
    n_top_genes=params.n_top_genes,
    min_mean=params.min_mean,
    max_mean=params.max_mean,
    min_disp=params.min_disp
)

# Store raw data and subset to highly variable genes
adata.raw = adata.copy()
adata = adata[:, adata.var.highly_variable]
print(f"  After HVG selection: {adata.shape}")

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

# Step 1: Feature auto-encoder
print("\n=== Step 1: Feature Auto-Encoder ===")
print(f"  Max epochs: {params.max_epoch_fae}")
print(f"  Learning rate: {params.lr_fae}")
print(f"  Hidden layers: {params.hidden_layers_fae}")

cd.encode_features(
    adata,
    D_encode_list=params.hidden_layers_fae + [params.latent_size_fae],
    D_decode_list=list(reversed(params.hidden_layers_fae)) + [adata.n_vars],
    max_epoch=params.max_epoch_fae,
    lr=params.lr_fae,
    device=params.device
)

print(f"Feature encoding complete! Saved to adata.obsm['X_fae']")

# Step 2: Create node_batch_mt (required by build_integration_loss_adj)
print("\n=== Step 2: Creating batch matrix ===")
if 'node_batch_mt' not in adata.obsm:
    adata.obsm['node_batch_mt'] = pd.get_dummies(adata.obs[params.batch_key]).to_numpy()
    print("Created node_batch_mt")
else:
    print("node_batch_mt already exists")

# Step 3: Build integration loss adjacency
print("\n=== Step 3: Building integration loss adjacency ===")
print(f"  k: {params.k}")

cd.inte.build_integration_loss_adj(
    adata,
    use_rep='X_fae',
    k=params.k,
    device=params.device
)

# Step 4: Build integration graph
print("\n=== Step 4: Building integration graph ===")
print(f"  n_edges_per_node: {params.n_edges_per_node}")
print(f"  k_mnn: {params.k_mnn}")

cd.inte.build_integration_graph(
    adata,
    batch_key=params.batch_key,
    use_rep='X_fae',
    n_edges_per_node=params.n_edges_per_node,
    k_mnn=params.k_mnn,
    device=params.device
)

print(f"Integration graph shape: {adata.uns['integration_edge_index'].shape}")

# Step 5: Run CellDiffusion integration
print("\n=== Step 5: Running CellDiffusion integration ===")
print(f"  Max epochs: {params.max_epoch_dif}")
print(f"  Learning rate: {params.lr_dif}")
print(f"  Num features: {params.num_features_diffusion}")
print(f"  Num heads: {params.num_heads_diffusion}")
print(f"  Num steps: {params.num_steps_diffusion}")
print(f"  Time increment: {params.time_increment_diffusion}")

cd.inte.integration_diffusion(
    adata,
    use_rep='X_fae',
    save_key='X_dif',
    max_epoch=params.max_epoch_dif,
    lr=params.lr_dif,
    num_features_diffusion=params.num_features_diffusion,
    num_heads_diffusion=params.num_heads_diffusion,
    num_steps_diffusion=params.num_steps_diffusion,
    time_increment_diffusion=params.time_increment_diffusion,
    device=params.device
)

print(f"CellDiffusion integration complete! Saved to adata.obsm['X_dif']")

# Step 6: Compute UMAP
print("\n=== Step 6: Computing UMAP ===")
print(f"  n_neighbors: {params.umap_n_neighbors}")
print(f"  min_dist: {params.umap_min_dist}")
print(f"  n_pcs: {params.umap_n_pcs}")

sc.pp.neighbors(adata, use_rep='X_dif', n_neighbors=params.umap_n_neighbors, n_pcs=params.umap_n_pcs)
sc.tl.umap(adata, min_dist=params.umap_min_dist)

# Rename to X_umap_dif for clarity
adata.obsm['X_umap_dif'] = adata.obsm['X_umap'].copy()
print(f"UMAP complete! Saved to adata.obsm['X_umap_dif']")

# Step 7: Leiden clustering
print("\n=== Step 7: Leiden clustering ===")
print(f"  Resolution: {params.leiden_resolution}")

sc.tl.leiden(adata, resolution=params.leiden_resolution, key_added='leiden_dif')
print(f"Leiden clustering complete! Saved to adata.obs['leiden_dif']")

# Save integrated data
print(f"\n=== Saving integrated data ===")
print(f"  Output: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

print("\n=== CellDiffusion pipeline complete! ===")
