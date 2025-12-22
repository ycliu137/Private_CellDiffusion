"""
Integration pipeline: build graph, add random edges, run diffusion, evaluate
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

# Load input data
input_h5ad = snakemake.input.h5ad
output_h5ad = snakemake.output.h5ad
metrics_csv = snakemake.output.metrics
params = snakemake.params

print(f"Loading encoded data from: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)

print(f"Data shape: {adata.shape}")

# Step 0: Create node_batch_mt (required by build_integration_loss_adj)
print("\n=== Creating node_batch_mt ===")
if 'node_batch_mt' not in adata.obsm:
    adata.obsm['node_batch_mt'] = pd.get_dummies(adata.obs[params.batch_key]).to_numpy()
    print("Created node_batch_mt")
else:
    print("node_batch_mt already exists")

# Step 1: Build integration loss adjacency
print("\n=== Building integration loss adjacency ===")
cd.inte.build_integration_loss_adj(
    adata,
    use_rep='X_fae',
    k=params.loss_adj_k,
    device=params.device
)

# Step 2: Build integration graph
print("\n=== Building integration graph ===")
cd.inte.build_integration_graph(
    adata,
    batch_key=params.batch_key,
    use_rep='X_fae',
    n_edges_per_node=params.n_edges_per_node,
    k_mnn=params.k_mnn,
    device=params.device
)

print(f"Integration graph shape: {adata.uns['integration_edge_index'].shape}")

# Step 3: Evaluate neighbor purity before adding random edges
print("\n=== Evaluating neighbor purity (before random edges) ===")
neighbor_purity_before = cd.eval.evaluate_neighbor_purity(
    adata,
    label_key=params.label_key,
    graph_key='integration_edge_index'
)
print(f"Neighbor purity (before): {neighbor_purity_before}")

# Step 4: Add random edges
print(f"\n=== Adding random edges (k_add={params.k_add}, seed={params.seed}) ===")
cd.eval.add_random_edges(
    adata,
    batch_key=params.batch_key,
    to_graph="integration_edge_index",
    k_add=params.k_add,
    seed=params.seed,
    device=params.device
)

print(f"Integration graph shape after adding random edges: {adata.uns['integration_edge_index'].shape}")

# Step 5: Run integration diffusion
print("\n=== Running integration diffusion ===")
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
    max_epoch=params.max_epoch,
    lr=params.lr,
    num_features_diffusion=params.num_features_diffusion,
    num_heads_diffusion=params.num_heads_diffusion,
    num_steps_diffusion=params.num_steps_diffusion,
    time_increment_diffusion=params.time_increment_diffusion,
    device=params.device
)

# Step 6: Evaluate MNN neighbor purity after integration
print("\n=== Evaluating MNN neighbor purity (after integration) ===")
neighbor_purity_after, _ = cd.eval.evaluate_mnn_neighbor_purity(
    adata,
    use_rep='X_dif',
    batch_key=params.batch_key,
    label_key=params.label_key,
    k_mnn=params.eval_k_mnn
)
print(f"Neighbor purity (after): {neighbor_purity_after}")

# Step 7: Compute UMAP for visualization
print("\n=== Computing UMAP ===")
sc.pp.neighbors(adata, use_rep='X_dif', n_neighbors=50, n_pcs=50)
sc.tl.umap(adata)

# Save integrated data
print(f"\n=== Saving integrated data ===")
print(f"Saving to: {output_h5ad}")
Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_h5ad)

# Save metrics
print(f"\n=== Saving metrics ===")
print(f"Saving metrics to: {metrics_csv}")
# Get k_value and k_add_value from params (passed from Snakefile)
k_value = getattr(params, 'k_value', params.loss_adj_k)
k_add_value = getattr(params, 'k_add_value', params.k_add)
# Write metrics with header (each combination writes to its own file)
with open(metrics_csv, mode="w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["loss_adj_k", "k_added", "neighbor_purity_before", "neighbor_purity_after"])
    writer.writerow([k_value, k_add_value, neighbor_purity_before, neighbor_purity_after])

print("\n=== Integration pipeline complete! ===")

