"""
SCIB evaluation: benchmark different k_add integration results using scib-metrics
"""
import sys
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless environments

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from scib_metrics.benchmark import Benchmarker
from scib_metrics.benchmark import BatchCorrection
from scib_metrics_new.benchmarked_plot_results_table import plot_results_table_new

# Get input and output from snakemake
input_h5ad = snakemake.input.h5ad
output_results = snakemake.output.results
output_table = snakemake.output.table
output_plot = snakemake.output.plot
params = snakemake.params

print(f"\n=== Loading aggregated adata ===")
print(f"Input file: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)
print(f"Data shape: {adata.shape}")

# Get all X_dif_kadd{k_add} keys from obsm
print(f"\n=== Extracting X_dif representations ===")
x_dif_keys = [key for key in adata.obsm.keys() if key.startswith('X_dif_kadd')]
x_dif_keys = sorted(x_dif_keys, key=lambda x: int(x.split('kadd')[1]))  # Sort by k_add value
print(f"Found {len(x_dif_keys)} X_dif representations:")
for key in x_dif_keys:
    print(f"  - {key}: shape {adata.obsm[key].shape}")

# Determine pre-integrated embedding key
# Use X_fae if available, otherwise use the first PCA representation
pre_integrated_key = None
if 'X_fae' in adata.obsm:
    pre_integrated_key = 'X_fae'
    print(f"\nUsing 'X_fae' as pre-integrated embedding")
elif 'X_pca' in adata.obsm:
    pre_integrated_key = 'X_pca'
    print(f"\nUsing 'X_pca' as pre-integrated embedding")
else:
    # Use the first available representation as fallback
    available_keys = [k for k in adata.obsm.keys() if not k.startswith('X_dif_kadd')]
    if available_keys:
        pre_integrated_key = available_keys[0]
        print(f"\nWarning: Using '{pre_integrated_key}' as pre-integrated embedding (fallback)")
    else:
        print("\nWarning: No suitable pre-integrated embedding found, will skip PCR comparison")

# Set up batch correction metrics
batch_corr = BatchCorrection(pcr_comparison=(pre_integrated_key is not None))

print(f"\n=== Running SCIB benchmark ===")
print(f"  Batch key: {params.batch_key}")
print(f"  Label key: {params.label_key}")
print(f"  Embedding keys: {x_dif_keys}")
if pre_integrated_key:
    print(f"  Pre-integrated key: {pre_integrated_key}")

# Create Benchmarker
bm = Benchmarker(
    adata,
    batch_key=params.batch_key,
    label_key=params.label_key,
    embedding_obsm_keys=x_dif_keys,
    pre_integrated_embedding_obsm_key=pre_integrated_key,
    batch_correction_metrics=batch_corr,
    n_jobs=params.n_jobs if hasattr(params, 'n_jobs') else 1,
)

# Run benchmark
print("\nStarting benchmark...")
bm.benchmark()

print("\n=== Benchmark complete ===")

# Save results
print(f"\n=== Saving results ===")

# Save results table as CSV
results_df = bm.get_results(min_max_scale=False)
print(f"Saving results table to: {output_table}")
Path(output_table).parent.mkdir(parents=True, exist_ok=True)
results_df.to_csv(output_table, index=True)
print(f"Results table shape: {results_df.shape}")
print("\nResults summary:")
print(results_df)

# Generate and save plots
print(f"\n=== Generating and saving plots ===")

# Assign custom plot function to benchmarker
bm.plot_results_table_new = plot_results_table_new

# Generate plot using custom plot function and save to file
print(f"Generating results table plot and saving to: {output_plot}")
Path(output_plot).parent.mkdir(parents=True, exist_ok=True)

try:
    fig = bm.plot_results_table_new(min_max_scale=False, show=False, save_fig=str(output_plot))
    if fig is not None:
        # Close the figure to free memory (plot_results_table_new already saved it)
        import matplotlib.pyplot as plt
        plt.close(fig)
    print("  Plot generated and saved successfully")
except Exception as e:
    print(f"  Warning: Could not generate plot with custom function: {e}")
    print("  Trying original plot function...")
    try:
        import matplotlib.pyplot as plt
        fig = bm.plot_results_table(min_max_scale=False)
        if fig is not None:
            plt.savefig(output_plot, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"  Original plot function executed and saved successfully")
        else:
            print(f"  Warning: Original plot function returned None")
    except Exception as e2:
        print(f"  Warning: Original plot function also failed: {e2}")
        print("  Continuing without plot...")

# Save full benchmarker object (pickle)
print(f"\nSaving benchmarker object to: {output_results}")
import pickle
Path(output_results).parent.mkdir(parents=True, exist_ok=True)
with open(output_results, 'wb') as f:
    pickle.dump(bm, f)

print("\n=== SCIB evaluation complete! ===")

