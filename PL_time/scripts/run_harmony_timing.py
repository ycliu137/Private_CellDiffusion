"""
Harmony integration with timing measurements.
Records dataset statistics, running time, and max CPU memory.
Harmony uses PyTorch; force CPU to avoid GPU conflicts when only one GPU is available.
"""
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # Force CPU (must be before torch/harmonypy import)

import sys
import json
import time
from pathlib import Path
from datetime import datetime

# Add project root and PL_time scripts to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(Path(__file__).parent))
from memory_monitor import start_cpu_monitor, stop_cpu_monitor

import scanpy as sc
import harmonypy
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
    "method": "Harmony",
    "timestamp": datetime.now().isoformat(),
    "steps": {}
}

# Initialize stats dictionary
stats_dict = {
    "dataset": Path(input_h5ad).parent.name,
    "method": "Harmony",
    "timestamp": datetime.now().isoformat()
}

print(f"\n=== Harmony Timing Pipeline ===")
print(f"Input: {input_h5ad}")
print(f"Output: {output_h5ad}")

# Start memory monitoring
_ = start_cpu_monitor()

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

adata = adata[:, adata.var.highly_variable]
print(f"Shape after HVG selection: {adata.shape}")

t_preprocess = time.time() - t0
timing_dict["steps"]["preprocessing"] = t_preprocess
print(f"Preprocessing time: {t_preprocess:.2f}s")

# ===== Step 1: Scaling =====
print(f"\n=== Step 1: Scaling ===")
t0 = time.time()

sc.pp.scale(adata, max_value=10)

t_scale = time.time() - t0
timing_dict["steps"]["scaling"] = t_scale
print(f"Scaling time: {t_scale:.2f}s")

# ===== Step 2: PCA =====
print(f"\n=== Step 2: PCA ===")
t0 = time.time()

sc.tl.pca(adata, svd_solver='arpack', n_comps=params.n_pcs)

t_pca = time.time() - t0
timing_dict["steps"]["pca"] = t_pca
print(f"PCA time: {t_pca:.2f}s")

# ===== Step 3: Harmony Integration =====
print(f"\n=== Step 3: Harmony Integration ===")
print(f"  Batch key: {params.batch_key}")
print(f"  Theta: {params.theta}")
print(f"  Lambda: {params.lambda_}")
print(f"  Sigma: {params.sigma}")

t0 = time.time()

# Run Harmony
ho = harmonypy.run_harmony(
    adata.obsm['X_pca'],
    adata.obs,
    vars_use=[params.batch_key],
    theta=params.theta,
    lamb=params.lambda_,
    sigma=params.sigma,
    nclust=params.nclust,
    max_iter_harmony=params.max_iter_harmony,
    verbose=False
)

# Update adata with Harmony embeddings
# ho.Z_corr is already (n_cells, n_pcs), no transpose needed
adata.obsm['X_harmony'] = ho.Z_corr
print(f"Harmony embedding shape: {adata.obsm['X_harmony'].shape}")

t_harmony = time.time() - t0
timing_dict["steps"]["harmony"] = t_harmony
print(f"Harmony integration time: {t_harmony:.2f}s")

# ===== Step 4: UMAP =====
print(f"\n=== Step 4: UMAP ===")
t0 = time.time()

sc.pp.neighbors(adata, use_rep='X_harmony', n_neighbors=15, n_pcs=params.n_pcs)
sc.tl.umap(adata, min_dist=0.3)
adata.obsm['X_umap_harmony'] = adata.obsm['X_umap'].copy()

t_umap = time.time() - t0
timing_dict["steps"]["umap"] = t_umap
print(f"UMAP time: {t_umap:.2f}s")

# ===== Step 5: Leiden =====
print(f"\n=== Step 5: Leiden Clustering ===")
t0 = time.time()

sc.tl.leiden(adata, resolution=params.leiden_resolution, key_added='leiden_harmony')

t_leiden = time.time() - t0
timing_dict["steps"]["leiden"] = t_leiden
print(f"Leiden time: {t_leiden:.2f}s")

# ===== Save results =====
print(f"\n=== Saving results ===")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)
print(f"Data saved to: {output_h5ad}")

# Calculate total time
def _step_duration(v):
    if isinstance(v, dict):
        return v.get("duration", 0)
    return v

total_time = sum(_step_duration(v) for v in timing_dict["steps"].values())
for key in timing_dict["steps"]:
    if isinstance(timing_dict["steps"][key], dict):
        timing_dict["steps"][key]["duration"] = timing_dict["steps"][key]["duration"] / 60
    else:
        timing_dict["steps"][key] = timing_dict["steps"][key] / 60
timing_dict["total_time"] = total_time / 60
print(f"\nTotal Harmony time: {total_time/60:.2f}min")

# Save timing results
Path(output_timing).parent.mkdir(parents=True, exist_ok=True)
with open(output_timing, 'w') as f:
    json.dump(timing_dict, f, indent=2)
print(f"Timing saved to: {output_timing}")

# Record max memory (Harmony is CPU-only)
stats_dict["max_cpu_memory_gb"] = stop_cpu_monitor()
stats_dict["max_gpu_memory_gb"] = None
print(f"Max CPU memory: {stats_dict['max_cpu_memory_gb']:.3f} GB" if stats_dict['max_cpu_memory_gb'] is not None else "Max CPU memory: N/A (psutil required)")

# Save statistics
Path(output_stats).parent.mkdir(parents=True, exist_ok=True)
with open(output_stats, 'w') as f:
    json.dump(stats_dict, f, indent=2)
print(f"Statistics saved to: {output_stats}")

print("\n=== Harmony timing pipeline complete! ===")
