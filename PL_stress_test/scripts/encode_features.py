"""
Feature encoding using autoencoder
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

print(f"Loading preprocessed data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

print(f"Data shape: {adata.shape}")

# CUDA synchronization: Wait for GPU to be available
if params.device == "cuda":
    print("\n=== Waiting for GPU availability ===")
    if torch.cuda.is_available():
        # Add small random delay to stagger simultaneous task starts
        initial_delay = random.uniform(0, 2)
        time.sleep(initial_delay)
        wait_interval = 0.5  # Check every 0.5 seconds
        waited = 0
        while True:
            try:
                # Try to create a small tensor on GPU to test availability
                test_tensor = torch.zeros(1, device='cuda')
                torch.cuda.synchronize()
                del test_tensor
                torch.cuda.empty_cache()
                print(f"GPU is available (waited {waited:.1f}s)")
                break
            except RuntimeError as e:
                if "busy" in str(e).lower() or "unavailable" in str(e).lower():
                    print(f"GPU is busy, waiting... ({waited:.1f}s)")
                    time.sleep(wait_interval)
                    waited += wait_interval
                else:
                    raise

# Encode features
print("Starting feature encoder...")
print(f"  Encoder dimensions: {params.D_encode_list}")
print(f"  Decoder dimensions: {params.D_decode_list}")
print(f"  Max epochs: {params.max_epoch}")
print(f"  Learning rate: {params.lr}")
print(f"  Device: {params.device}")

cd.encode_features(
    adata,
    D_encode_list=params.D_encode_list,
    D_decode_list=params.D_decode_list,
    max_epoch=params.max_epoch,
    lr=params.lr,
    device=params.device
)

print("Feature encoding complete!")

# Save encoded data
print(f"Saving encoded data to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

print("Encoding step complete!")

