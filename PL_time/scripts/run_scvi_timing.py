"""
scVI integration with timing measurements.
Records dataset statistics and running time.
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
import scvi
import torch
import pandas as pd

# Suppress scvi logging
scvi.settings.verbosity = 0

# Load input data
input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
output_timing = snakemake.output.timing
output_stats = snakemake.output.stats
params = snakemake.params

# Initialize timing dictionary
timing_dict = {
    "dataset": Path(input_h5ad).parent.name,
    "method": "scVI",
    "timestamp": datetime.now().isoformat(),
    "steps": {}
}

# Initialize stats dictionary
stats_dict = {
    "dataset": Path(input_h5ad).parent.name,
    "method": "scVI",
    "timestamp": datetime.now().isoformat()
}

print(f"\n=== scVI Timing Pipeline ===")
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
    # scVI requires raw counts (not log-normalized)
    # Do NOT normalize or log1p
else:
    print("Data already in raw counts, using as-is for scVI")

print(f"Finding highly variable genes (batch-aware): n_top_genes={params.n_top_genes}")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=params.n_top_genes,
    min_mean=params.min_mean,
    max_mean=params.max_mean,
    min_disp=params.min_disp,
    batch_key=params.batch_key,
    flavor="seurat_v3"
)

adata = adata[:, adata.var.highly_variable]
adata = adata.copy()  # Convert view to copy for scVI compatibility
print(f"Shape after HVG selection: {adata.shape}")

t_preprocess = time.time() - t0
timing_dict["steps"]["preprocessing"] = t_preprocess
print(f"Preprocessing time: {t_preprocess:.2f}s")

# ===== Step 1: Setup scVI =====
print(f"\n=== Step 1: Setup scVI ===")
t0 = time.time()

# Initialize scVI model
scvi.model.SCVI.setup_anndata(adata, batch_key=params.batch_key)

t_setup = time.time() - t0
timing_dict["steps"]["setup_scvi"] = t_setup
print(f"Setup time: {t_setup:.2f}s")

# ===== Step 2: Train scVI =====
print(f"\n=== Step 2: Training scVI ===")
print(f"  n_layers: {params.n_layers}")
print(f"  n_latent: {params.n_latent}")
print(f"  max_epochs: {params.max_epochs}")
print(f"  early_stopping: {params.early_stopping}")

t0 = time.time()

model = scvi.model.SCVI(
    adata,
    n_layers=params.n_layers,
    n_latent=params.n_latent,
    dispersion="gene"
)

# Prepare training kwargs
train_kwargs = {
    'max_epochs': params.max_epochs,
    'early_stopping': params.early_stopping,
    'accelerator': 'cpu',
    'progress_bar': False
}

if params.early_stopping:
    train_kwargs['early_stopping_patience'] = params.early_stopping_patience

model.train(**train_kwargs)

# Get latent representation
latent = model.get_latent_representation()
adata.obsm['X_scvi'] = latent
print(f"scVI embedding shape: {adata.obsm['X_scvi'].shape}")

t_train = time.time() - t0
timing_dict["steps"]["train_scvi"] = t_train
print(f"Training time: {t_train:.2f}s")

# Note: scVI latent representation is already low-dimensional (n_latent=10)
# No additional PCA needed

# ===== Step 3: UMAP (on scVI latent) =====
print(f"\n=== Step 3: UMAP (on scVI latent) ===")
t0 = time.time()

sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=15)
sc.tl.umap(adata, min_dist=0.3)
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

t_umap = time.time() - t0
timing_dict["steps"]["umap"] = t_umap
print(f"UMAP time: {t_umap:.2f}s")

# ===== Step 4: Leiden Clustering ===
print(f"\n=== Step 4: Leiden Clustering ===")
t0 = time.time()

sc.tl.leiden(adata, resolution=params.leiden_resolution, key_added='leiden_scvi')

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
for key in timing_dict["steps"]:
    timing_dict["steps"][key] = timing_dict["steps"][key] / 60
timing_dict["total_time"] = total_time / 60
print(f"\nTotal scVI time: {total_time/60:.2f}min")

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

print("\n=== scVI timing pipeline complete! ===")
