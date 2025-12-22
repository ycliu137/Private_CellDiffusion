"""
Feature encoding using autoencoder
"""
import sys
from pathlib import Path
import scanpy as sc
import celldiffusion as cd

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Load input data
input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
params = snakemake.params

print(f"Loading preprocessed data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

print(f"Data shape: {adata.shape}")

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

