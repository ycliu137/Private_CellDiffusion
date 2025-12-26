######## snakemake preamble start (automatically inserted, do not edit) ########
import sys;sys.path.extend(['/projectnb/czproj/packages/ycliu/envs/dif_snake_env/lib/python3.12/site-packages', '/projectnb/czproj/Algorithms/ycliu/Private_CellDiffusion/PL_stress_test', '/projectnb/czproj/packages/ycliu/envs/dif_snake_env/bin', '/projectnb/czproj/packages/ycliu/envs/dif_snake_env/lib/python3.12', '/projectnb/czproj/packages/ycliu/envs/dif_snake_env/lib/python3.12/lib-dynload', '/projectnb/czproj/packages/ycliu/envs/dif_snake_env/lib/python3.12/site-packages', '/usr2/postdoc/ycliu137/.cache/snakemake/snakemake/source-cache/runtime-cache/tmpz20pelbm/file/projectnb/czproj/Algorithms/ycliu/Private_CellDiffusion/PL_stress_test/scripts', '/projectnb/czproj/Algorithms/ycliu/Private_CellDiffusion/PL_stress_test/scripts']);import pickle;from snakemake import script;script.snakemake = pickle.loads(b"\x80\x04\x95,\n\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8cg/projectnb/czproj/Algorithms/ycliu/Private_CellDiffusion/data/outputs/stress_tests_pbmc10k/encoded.h5ad\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x04h5ad\x94K\x00N\x86\x94s\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x12h\x06\x8c\x0eAttributeGuard\x94\x93\x94)\x81\x94}\x94\x8c\x04name\x94h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbh\x0eh\nub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c\x82/projectnb/czproj/Algorithms/ycliu/Private_CellDiffusion/data/outputs/stress_tests_pbmc10k/scDiffusion_integration_k50_kadd70.h5ad\x94\x8cu/projectnb/czproj/Algorithms/ycliu/Private_CellDiffusion/data/outputs/stress_tests_pbmc10k/metrics_log_k50_kadd70.csv\x94e}\x94(h\x0c}\x94(h\x0eK\x00N\x86\x94\x8c\x07metrics\x94K\x01N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x15)\x81\x94}\x94h\x18h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbh\x0eh\x1fh$h ub\x8c\r_params_store\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(K2K2K2KFK*M\xd0\x07G?PbM\xd2\xf1\xa9\xfcK2K\x08K\x08G?\xc9\x99\x99\x99\x99\x99\x9a\x8c\x06labels\x94\x8c\x05batch\x94K2\x8c\x04cuda\x94K2KFe}\x94(h\x0c}\x94(\x8c\x10n_edges_per_node\x94K\x00N\x86\x94\x8c\x05k_mnn\x94K\x01N\x86\x94\x8c\nloss_adj_k\x94K\x02N\x86\x94\x8c\x05k_add\x94K\x03N\x86\x94\x8c\x04seed\x94K\x04N\x86\x94\x8c\tmax_epoch\x94K\x05N\x86\x94\x8c\x02lr\x94K\x06N\x86\x94\x8c\x16num_features_diffusion\x94K\x07N\x86\x94\x8c\x13num_heads_diffusion\x94K\x08N\x86\x94\x8c\x13num_steps_diffusion\x94K\tN\x86\x94\x8c\x18time_increment_diffusion\x94K\nN\x86\x94\x8c\tlabel_key\x94K\x0bN\x86\x94\x8c\tbatch_key\x94K\x0cN\x86\x94\x8c\neval_k_mnn\x94K\rN\x86\x94\x8c\x06device\x94K\x0eN\x86\x94\x8c\x07k_value\x94K\x0fN\x86\x94\x8c\x0bk_add_value\x94K\x10N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x15)\x81\x94}\x94h\x18h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbh4K2h6K2h8K2h:KFh<K*h>M\xd0\x07h@G?PbM\xd2\xf1\xa9\xfchBK2hDK\x08hFK\x08hHG?\xc9\x99\x99\x99\x99\x99\x9ahJh/hLh0hNK2hPh1hRK2hTKFub\x8c\r_params_types\x94}\x94\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94(\x8c\x0250\x94\x8c\x0270\x94e}\x94(h\x0c}\x94(\x8c\x01k\x94K\x00N\x86\x94\x8c\x05k_add\x94K\x01N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x15)\x81\x94}\x94h\x18h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbhehah:hbub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c\x1d/scratch/2460558.1.czcb-buyin\x94K\x01e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94\x8c\x03gpu\x94K\x03N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x15)\x81\x94}\x94h\x18h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbhvK\x01hxK\x01hzhsh|K\x01ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x15)\x81\x94}\x94h\x18h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbub\x8c\x06config\x94}\x94(\x8c\x04data\x94}\x94(\x8c\ninput_path\x94\x8cT/projectnb/czlab/A00/ycliu/scRNA_integration_benchmark_datasets/PBMC10K/pbmc10k.h5ad\x94\x8c\noutput_dir\x94\x8c$../data/outputs/stress_tests_pbmc10k\x94u\x8c\npreprocess\x94}\x94(\x8c\tmin_cells\x94K\x03\x8c\ntarget_sum\x94M\x10'\x8c\x0bn_top_genes\x94M\xd0\x07\x8c\x08min_mean\x94G?\x89\x99\x99\x99\x99\x99\x9a\x8c\x08max_mean\x94K\t\x8c\x08min_disp\x94G?\xe0\x00\x00\x00\x00\x00\x00u\x8c\x0ffeature_encoder\x94}\x94(\x8c\rD_encode_list\x94]\x94(M\xd0\x07M,\x01K2e\x8c\rD_decode_list\x94]\x94(K2M,\x01M\xd0\x07e\x8c\tmax_epoch\x94M\xd0\x07\x8c\x02lr\x94G?PbM\xd2\xf1\xa9\xfcu\x8c\x11integration_graph\x94}\x94(\x8c\x10n_edges_per_node\x94K2\x8c\x05k_mnn\x94K2u\x8c\x14integration_loss_adj\x94}\x94heK2s\x8c\x0crandom_edges\x94}\x94(\x8c\x05k_add\x94]\x94(K\x00K\nK\x14K\x1eK(K2K<KFKPKZKde\x8c\x04seed\x94K*u\x8c\x15integration_diffusion\x94}\x94(\x8c\tmax_epoch\x94M\xd0\x07\x8c\x02lr\x94G?PbM\xd2\xf1\xa9\xfc\x8c\x16num_features_diffusion\x94K2\x8c\x13num_heads_diffusion\x94K\x08\x8c\x13num_steps_diffusion\x94K\x08\x8c\x18time_increment_diffusion\x94G?\xc9\x99\x99\x99\x99\x99\x9au\x8c\nevaluation\x94}\x94(\x8c\tlabel_key\x94h/\x8c\tbatch_key\x94h0\x8c\x05k_mnn\x94K2u\x8c\x04umap\x94}\x94(\x8c\x0bn_neighbors\x94K2\x8c\x05n_pcs\x94K2u\x8c\x06device\x94h1u\x8c\x04rule\x94\x8c\tintegrate\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8cO/projectnb/czproj/Algorithms/ycliu/Private_CellDiffusion/PL_stress_test/scripts\x94ub.");del script;from snakemake.logging import logger;from snakemake.script import snakemake;__real_file__ = __file__; __file__ = '/projectnb/czproj/Algorithms/ycliu/Private_CellDiffusion/PL_stress_test/scripts/integrate.py';
######## snakemake preamble end #########
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

