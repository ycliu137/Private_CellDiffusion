"""
Run GCN integration with varying num_layers_gcn
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
print(f"Number of GCN layers: {params.num_layers_gcn}")

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

# Run GCN integration
print("\n=== Running GCN integration ===")
print(f"  Max epochs: {params.max_epoch}")
print(f"  Learning rate: {params.lr}")
print(f"  Num features GCN: {params.num_features_gcn}")
print(f"  Num layers GCN: {params.num_layers_gcn}")
print(f"  Dropout: {params.dropout}")
print(f"  Device: {params.device}")

cd.inte.integration_gcn(
    adata,
    use_rep='X_fae',
    save_key='X_gcn',  # Save as X_gcn
    max_epoch=params.max_epoch,
    lr=params.lr,
    num_features_gcn=params.num_features_gcn,
    num_layers_gcn=params.num_layers_gcn,
    dropout=params.dropout,
    encoder=None,  # Use default
    decoder=None,  # Use default
    save_model=False,
    load_model_state=False,
    loss_reduction="sum",
    device=params.device
)

print(f"\nGCN integration complete! Embeddings saved in adata.obsm['X_gcn']")

# Save integrated data
print(f"\n=== Saving integrated data ===")
print(f"Saving to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

print("\n=== GCN integration complete! ===")

