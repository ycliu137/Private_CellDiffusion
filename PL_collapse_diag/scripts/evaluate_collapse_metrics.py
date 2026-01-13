"""
Evaluate collapse diagnostic metrics: intrinsic_dimension_knn, intrinsic_dimension_at_variance_percentage, variance_explained_by_embedding
"""
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import pandas as pd
from sc_evaluation.collapse_diag import (
    intrinsic_dimension_knn,
    intrinsic_dimension_at_variance_percentage,
    variance_explained_by_embedding
)

# Get input/output from snakemake
input_h5ad = snakemake.input.h5ad
output_csv = snakemake.output.csv
params = snakemake.params

print(f"Loading integrated data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

print(f"Data shape: {adata.shape}")
print(f"Method: {params.method}")
print(f"Network layers: {params.network_layers}")

# Determine which representation to use
if params.method == "CellDiffusion":
    use_rep = 'X_dif'
    print(f"Using representation: {use_rep}")
elif params.method == "GCN":
    use_rep = 'X_gcn'
    print(f"Using representation: {use_rep}")
else:
    raise ValueError(f"Unknown method: {params.method}")

# Check if representation exists
if use_rep not in adata.obsm:
    raise ValueError(f"Representation '{use_rep}' not found in adata.obsm. Available keys: {list(adata.obsm.keys())}")

# Evaluate metrics
print(f"\n=== Evaluating collapse diagnostic metrics ===")
print(f"  Use rep: {use_rep}")
print(f"  k for intrinsic_dimension_knn: {params.k}")
print(f"  variance_percentage: {params.variance_percentage}")

# 1. Intrinsic dimension (KNN)
print(f"\n  Computing intrinsic_dimension_knn...")
intrinsic_dim_knn = intrinsic_dimension_knn(adata, use_rep=use_rep, k=params.k)
print(f"    Intrinsic dimension (KNN): {intrinsic_dim_knn:.4f}")

# 2. Intrinsic dimension (variance percentage)
print(f"\n  Computing intrinsic_dimension_at_variance_percentage...")
intrinsic_dim_var = intrinsic_dimension_at_variance_percentage(
    adata, 
    use_rep=use_rep, 
    variance_percentage=params.variance_percentage
)
print(f"    Intrinsic dimension (variance {params.variance_percentage*100}%): {intrinsic_dim_var:.2f}")

# 3. Variance explained by embedding
# Check if adata.X exists (required for variance_explained_by_embedding)
variance_explained = None
if adata.X is not None:
    print(f"\n  Computing variance_explained_by_embedding...")
    try:
        variance_explained = variance_explained_by_embedding(adata, use_rep=use_rep)
        print(f"    Variance explained: {variance_explained:.4f}")
    except Exception as e:
        print(f"    Warning: Failed to compute variance_explained_by_embedding: {e}")
        variance_explained = None
else:
    print(f"\n  Warning: adata.X is None, skipping variance_explained_by_embedding")

# Save results to CSV
print(f"\n=== Saving results ===")
results = pd.DataFrame({
    'method': [params.method],
    'network_layers': [params.network_layers],
    'intrinsic_dimension_knn': [intrinsic_dim_knn],
    'intrinsic_dimension_variance': [intrinsic_dim_var],
    'variance_explained': [variance_explained],
    'k': [params.k],
    'variance_percentage': [params.variance_percentage]
})

print(f"Saving to: {output_csv}")
Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
results.to_csv(output_csv, index=False)

print("Evaluation complete!")

