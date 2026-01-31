"""
CellDiffusion integration with timing measurements.
Records dataset statistics and running time for each step.
"""
import sys
import json
import time
from pathlib import Path
from datetime import datetime

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import celldiffusion as cd
import torch
import random
import pandas as pd

# Load input data
input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
output_timing = snakemake.output.timing
output_stats = snakemake.output.stats
params = snakemake.params

# Initialize timing dictionary
timing_dict = {
    "dataset": Path(input_h5ad).parent.name,
    "method": "CellDiffusion",
    "timestamp": datetime.now().isoformat(),
    "steps": {}
}

# Initialize stats dictionary
stats_dict = {
    "dataset": Path(input_h5ad).parent.name,
    "method": "CellDiffusion",
    "timestamp": datetime.now().isoformat()
}

print(f"\n=== CellDiffusion Timing Pipeline ===")
print(f"Input: {input_h5ad}")
print(f"Output: {output_h5ad}")

# Load data
t0 = time.time()
adata = sc.read_h5ad(input_h5ad)
print(f"Data shape: {adata.shape}")

# Record dataset statistics
stats_dict["n_cells"] = adata.n_obs
stats_dict["n_genes"] = adata.n_vars
if params.batch_key in adata.obs.columns:
    stats_dict["n_batches"] = adata.obs[params.batch_key].nunique()
    print(f"Batches: {stats_dict['n_batches']}")
else:
    stats_dict["n_batches"] = 1
    print("No batch key found, assuming single batch")

timing_dict["steps"]["load_data"] = {"start": t0, "duration": time.time() - t0}

# ===== Step 0: Preprocessing =====
print(f"\n=== Step 0: Preprocessing ===")
t0 = time.time()

if not params.normalized_data:
    print(f"Filtering genes with min_cells={params.min_cells}")
    sc.pp.filter_genes(adata, min_cells=params.min_cells)
    
    print(f"Normalizing with target_sum={params.target_sum}")
    sc.pp.normalize_total(adata, target_sum=params.target_sum)
    sc.pp.log1p(adata)
else:
    print("Data already normalized, skipping preprocessing")

print(f"Finding highly variable genes: n_top_genes={params.n_top_genes}")
sc.pp.highly_variable_genes(
    adata, 
    n_top_genes=params.n_top_genes,
    min_mean=params.min_mean,
    max_mean=params.max_mean,
    min_disp=params.min_disp
)

adata.raw = adata.copy()
adata = adata[:, adata.var.highly_variable]
print(f"Shape after preprocessing: {adata.shape}")

t_preprocess = time.time() - t0
timing_dict["steps"]["preprocessing"] = t_preprocess
print(f"Preprocessing time: {t_preprocess:.2f}s")

# ===== Step 1: PCA =====
print(f"\n=== Step 1: PCA ===")
t0 = time.time()
sc.tl.pca(adata, svd_solver='arpack', n_comps=params.pca_n_comps)
t_pca = time.time() - t0
timing_dict["steps"]["pca"] = t_pca
print(f"PCA time: {t_pca:.2f}s")

# ===== Step 2: Feature Auto-Encoder =====
print(f"\n=== Step 2: Feature Auto-Encoder (FAE) ===")
t0 = time.time()

D_encode_list = [adata.n_vars] + params.hidden_layers_fae + [params.latent_size_fae]
D_decode_list = [params.latent_size_fae] + list(reversed(params.hidden_layers_fae)) + [adata.n_vars]

print(f"  Encoder architecture: {D_encode_list}")
print(f"  Decoder architecture: {D_decode_list}")

cd.encode_features(
    adata,
    D_encode_list=D_encode_list,
    D_decode_list=D_decode_list,
    max_epoch=params.max_epoch_fae,
    lr=params.lr_fae,
    device=params.device,
    batch_size=32
)

t_fae = time.time() - t0
timing_dict["steps"]["feature_autoencoder"] = t_fae
print(f"Feature Auto-Encoder time: {t_fae:.2f}s")

# ===== Step 3: Graph Building =====
print(f"\n=== Step 3: Graph Building ===")
t0 = time.time()

adata.obsm['X_neighbors'] = adata.obsm['X_fae'].copy()

cd.build_graph(
    adata,
    use_rep='X_fae',
    k=params.k,
    n_edges_per_node=params.n_edges_per_node,
    k_mnn=params.k_mnn,
    seed=0
)

t_graph = time.time() - t0
timing_dict["steps"]["graph_building"] = t_graph
print(f"Graph building time: {t_graph:.2f}s")

# ===== Step 4: Diffusion =====
print(f"\n=== Step 4: Graph Diffusion ===")
t0 = time.time()

cd.graph_diffusion(
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

t_diffusion = time.time() - t0
timing_dict["steps"]["graph_diffusion"] = t_diffusion
print(f"Graph diffusion time: {t_diffusion:.2f}s")

# ===== Step 5: UMAP =====
print(f"\n=== Step 5: UMAP ===")
t0 = time.time()

sc.pp.neighbors(adata, use_rep='X_dif', n_neighbors=params.umap_n_neighbors, n_pcs=params.umap_n_pcs)
sc.tl.umap(adata, min_dist=params.umap_min_dist)
adata.obsm['X_umap_dif'] = adata.obsm['X_umap'].copy()

t_umap = time.time() - t0
timing_dict["steps"]["umap"] = t_umap
print(f"UMAP time: {t_umap:.2f}s")

# ===== Step 6: Leiden =====
print(f"\n=== Step 6: Leiden Clustering ===")
t0 = time.time()

sc.tl.leiden(adata, resolution=params.leiden_resolution, key_added='leiden_dif')

t_leiden = time.time() - t0
timing_dict["steps"]["leiden"] = t_leiden
print(f"Leiden time: {t_leiden:.2f}s")

# ===== Save results =====
print(f"\n=== Saving results ===")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)
print(f"Data saved to: {output_h5ad}")

# Calculate total time
total_time = sum(timing_dict["steps"].values())
timing_dict["total_time"] = total_time
print(f"\nTotal CellDiffusion time: {total_time:.2f}s")

# Save timing results
Path(output_timing).parent.mkdir(parents=True, exist_ok=True)
with open(output_timing, 'w') as f:
    json.dump(timing_dict, f, indent=2)
print(f"Timing saved to: {output_timing}")

# Save statistics
Path(output_stats).parent.mkdir(parents=True, exist_ok=True)
with open(output_stats, 'w') as f:
    json.dump(stats_dict, f, indent=2)
print(f"Statistics saved to: {output_stats}")

print("\n=== CellDiffusion timing pipeline complete! ===")
