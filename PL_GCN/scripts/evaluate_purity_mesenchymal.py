"""
Evaluate neighbor purity for Mesenchymal lineage cells only using evaluate_knn_neighbor_purity
"""
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import celldiffusion as cd
import pandas as pd

# Get input/output from snakemake
input_h5ad = snakemake.input.h5ad
output_csv = snakemake.output.csv
params = snakemake.params

print(f"Loading integrated data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

print(f"Original data shape: {adata.shape}")
print(f"Method: {params.method}")

# Filter for Mesenchymal lineage cells
print(f"\n=== Filtering for Mesenchymal lineage cells ===")
lineage_key = 'lineage'
if lineage_key not in adata.obs.columns:
    raise ValueError(f"Column '{lineage_key}' not found in adata.obs. Available columns: {list(adata.obs.columns)}")

original_n_cells = adata.shape[0]
adata = adata[adata.obs[lineage_key] == 'Mesenchymal'].copy()
print(f"Filtered data shape: {adata.shape}")
print(f"Number of Mesenchymal cells: {adata.shape[0]} (out of {original_n_cells} total cells)")

if adata.shape[0] == 0:
    raise ValueError("No Mesenchymal lineage cells found in the data")

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

# Evaluate neighbor purity
print(f"\n=== Evaluating KNN neighbor purity (Mesenchymal lineage only) ===")
print(f"  Use rep: {use_rep}")
print(f"  Label key: {params.label_key}")
print(f"  k: {params.k} (number of nearest neighbors for KNN graph)")

purity, edge_index = cd.eval.evaluate_knn_neighbor_purity(
    adata,
    use_rep=use_rep,
    label_key=params.label_key,
    k=params.k  # k: number of nearest neighbors for KNN graph construction
)

print(f"  Neighbor purity: {purity:.4f}")

# Save results to CSV
print(f"\n=== Saving results ===")
results = pd.DataFrame({
    'method': [params.method],
    'neighbor_purity': [purity],
    'k': [params.k]
})

print(f"Saving to: {output_csv}")
Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
results.to_csv(output_csv, index=False)

print("Evaluation complete!")

