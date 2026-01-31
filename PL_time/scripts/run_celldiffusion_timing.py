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
    device=params.device
)

t_fae = time.time() - t0
timing_dict["steps"]["feature_autoencoder"] = t_fae
print(f"Feature Auto-Encoder time: {t_fae:.2f}s")

# ===== Step 3: Graph Building =====
print(f"\n=== Step 3: Graph Building ===")
t0 = time.time()

# Create node_batch_mt (required by build_integration_loss_adj)
if 'node_batch_mt' not in adata.obsm:
    adata.obsm['node_batch_mt'] = pd.get_dummies(adata.obs[params.batch_key]).to_numpy()
    print("Created node_batch_mt")

# Build integration loss adjacency
cd.inte.build_integration_loss_adj(
    adata,
    use_rep='X_fae',
    k=params.k,
    device=params.device
)

# Build integration graph
cd.inte.build_integration_graph(
    adata,
    batch_key=params.batch_key,
    use_rep='X_fae',
    n_edges_per_node=params.n_edges_per_node,
    k_mnn=params.k_mnn,
    device=params.device
)

t_graph = time.time() - t0
timing_dict["steps"]["graph_building"] = t_graph
print(f"Graph building time: {t_graph:.2f}s")

# ===== Step 4: Diffusion =====
print(f"\n=== Step 4: Graph Diffusion ===")
t0 = time.time()

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

t_diffusion = time.time() - t0
timing_dict["steps"]["graph_diffusion"] = t_diffusion
print(f"Graph diffusion time: {t_diffusion:.2f}s")

# ===== Save results =====
print(f"\n=== Saving results ===")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)
print(f"Data saved to: {output_h5ad}")

 # Calculate total time (handle nested step dicts like {"duration": ...})
def _step_duration(v):
    if isinstance(v, dict):
        return v.get("duration", 0)
    return v

total_time = t_fae + t_graph + t_diffusion
timing_dict["total_time"] = total_time
print(f"\nTotal CellDiffusion time: {total_time/60:.2f}min")

# Convert all timings to minutes
for key in timing_dict["steps"]:
    if isinstance(timing_dict["steps"][key], dict):
        timing_dict["steps"][key]["duration"] = timing_dict["steps"][key]["duration"] / 60
    else:
        timing_dict["steps"][key] = timing_dict["steps"][key] / 60
timing_dict["total_time"] = timing_dict["total_time"] / 60

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
