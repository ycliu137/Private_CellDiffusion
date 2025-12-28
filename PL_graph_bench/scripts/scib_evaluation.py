"""
SCIB evaluation: benchmark different graph building method integration results
"""
import sys
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import pickle
import traceback
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless environments

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from scib_metrics.benchmark import Benchmarker
from scib_metrics.benchmark import BatchCorrection

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

# Get all X_dif_{method} keys from obsm
print(f"\n=== Extracting X_dif representations ===")
x_dif_keys = [key for key in adata.obsm.keys() if key.startswith('X_dif_')]
x_dif_keys = sorted(x_dif_keys)  # Sort alphabetically
print(f"Found {len(x_dif_keys)} X_dif representations:")
for key in x_dif_keys:
    print(f"  - {key}: shape {adata.obsm[key].shape}")

if len(x_dif_keys) == 0:
    raise ValueError("No X_dif representations found in adata.obsm")

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
    available_keys = [k for k in adata.obsm.keys() if not k.startswith('X_dif_')]
    if available_keys:
        pre_integrated_key = available_keys[0]
        print(f"\nWarning: Using '{pre_integrated_key}' as pre-integrated embedding (fallback)")
    else:
        print("\nWarning: No suitable pre-integrated embedding found, will skip PCR comparison")

# Set up batch correction metrics
batch_corr = BatchCorrection(pcr_comparison=False)

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

plot_saved = False
try:
    # Try to use custom plot function if available
    try:
        from scib_metrics_new.benchmarked_plot_results_table import plot_results_table_new
        bm.plot_results_table_new = plot_results_table_new
        fig = bm.plot_results_table_new(min_max_scale=False, show=False, save_fig=str(output_plot))
        if fig is not None:
            import matplotlib.pyplot as plt
            plt.close(fig)
        if Path(output_plot).exists():
            print("  Plot generated and saved successfully using custom function")
            plot_saved = True
        else:
            print("  Warning: Custom function executed but file not found, trying alternative...")
    except ImportError:
        print("  Custom plot function not available, trying original function...")
    except Exception as e:
        print(f"  Warning: Custom plot function failed: {e}")
        print("  Trying original plot function...")
    
    # Fallback to original plot function if custom failed
    if not plot_saved:
        import matplotlib.pyplot as plt
        try:
            fig = bm.plot_results_table(min_max_scale=False)
            if fig is not None:
                plt.savefig(output_plot, dpi=300, bbox_inches='tight')
                plt.close(fig)
                if Path(output_plot).exists():
                    print("  Plot generated and saved successfully using original function")
                    plot_saved = True
                else:
                    print("  Warning: Original function executed but file not found")
            else:
                print("  Warning: Plot function did not return a figure")
        except Exception as e:
            print(f"  Error in original plot function: {e}")
            traceback.print_exc()
            
except Exception as e:
    print(f"  Error: Could not generate plot: {e}")
    traceback.print_exc()

# Ensure output file exists (create empty file if plotting failed)
if not Path(output_plot).exists():
    print(f"  Creating empty plot file as placeholder...")
    Path(output_plot).parent.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, 'Plot generation failed.\nCheck logs for details.', 
                ha='center', va='center', fontsize=12, 
                transform=ax.transAxes)
        ax.set_title("Plot Generation Failed")
        ax.axis('off')
        plt.savefig(output_plot, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Placeholder plot saved to {output_plot}")
    except Exception as e:
        print(f"  Critical Error: Could not create placeholder plot: {e}")
        # If even placeholder fails, create a minimal empty file
        with open(output_plot, 'w') as f:
            f.write("Plot generation failed. Check logs.")

# Save full benchmarker object (pickle)
print(f"\nSaving benchmarker object to: {output_results}")
Path(output_results).parent.mkdir(parents=True, exist_ok=True)
with open(output_results, 'wb') as f:
    pickle.dump(bm, f)

print("\n=== SCIB evaluation complete! ===")

