"""
Run CellDiffusion integration with specific n_edges_per_node
"""
import sys
from pathlib import Path

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

print(f"Loading data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

print(f"Data shape: {adata.shape}")
print(f"Integration for n_edges_per_node={params.n_edges_per_node}")

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

# Run CellDiffusion integration
print("\n=== Running CellDiffusion integration ===")
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

print(f"\nCellDiffusion integration complete! Embeddings saved in adata.obsm['X_dif']")

# Save integrated data
print(f"\n=== Saving integrated data ===")
print(f"Saving to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

print(f"\n=== CellDiffusion integration complete for n_edges_per_node={params.n_edges_per_node}! ===")

