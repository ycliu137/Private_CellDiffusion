"""
Evaluate neighbor purity on integration graphs built from X_fae and X_dif.
"""
import sys
from pathlib import Path
import csv
import pandas as pd

# Add project root to path first, before importing celldiffusion
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import celldiffusion as cd

# Load input and output paths
input_h5ad = snakemake.input.h5ad
output_csv = snakemake.output.csv
params = snakemake.params

print(f"Loading data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

# Resolve dataset name
dataset_name = getattr(snakemake.wildcards, "dataset", Path(input_h5ad).parent.name)

# Validate required embeddings
for rep in ("X_fae", "X_dif"):
    if rep not in adata.obsm:
        raise ValueError(f"Missing required embedding {rep} in adata.obsm")


def build_graph_and_purity(use_rep: str, graph_key: str) -> float:
    print(f"\n=== Building integration graph for {use_rep} ===")
    cd.inte.build_integration_graph(
        adata,
        batch_key=params.batch_key,
        use_rep=use_rep,
        n_edges_per_node=params.n_edges_per_node,
        k_mnn=params.k_mnn,
        device=params.device,
    )
    adata.uns[graph_key] = adata.uns["integration_edge_index"].copy()

    purity = cd.eval.evaluate_neighbor_purity(
        adata,
        label_key=params.label_key,
        graph_key=graph_key,
    )
    print(f"Neighbor purity ({use_rep}): {purity}")
    return purity


purity_x_fae = build_graph_and_purity("X_fae", "integration_edge_index_X_fae")
purity_x_dif = build_graph_and_purity("X_dif", "integration_edge_index_X_dif")

# Save metrics
print(f"\n=== Saving results ===")
print(f"Saving to: {output_csv}")
Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
with open(output_csv, mode="w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["dataset", "neighbor_purity_X_fae", "neighbor_purity_X_dif"])
    writer.writerow([dataset_name, purity_x_fae, purity_x_dif])

print("Done.")
