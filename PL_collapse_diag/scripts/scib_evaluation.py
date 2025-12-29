"""
SCIB evaluation: benchmark all CellDiffusion and GCN embeddings across different network layers
"""
import sys
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import pickle
import traceback
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for plotting

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from scib_metrics.benchmark import Benchmarker
from scib_metrics.benchmark import BatchCorrection

# Get input and output from snakemake
try:
    input_h5ad = snakemake.input.h5ad
    output_results = snakemake.output.results
    output_table = snakemake.output.table
    output_plot = snakemake.output.plot
    params = snakemake.params
except Exception as e:
    print(f"Error getting snakemake inputs/outputs: {e}")
    traceback.print_exc()
    sys.exit(1)

try:
    print(f"\n=== Loading aggregated adata ===")
    print(f"Input file: {input_h5ad}")
    if not Path(input_h5ad).exists():
        raise FileNotFoundError(f"Input file not found: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    print(f"Data shape: {adata.shape}")
    print(f"Available obsm keys: {list(adata.obsm.keys())}")
except Exception as e:
    print(f"Error loading adata: {e}")
    traceback.print_exc()
    sys.exit(1)

# Get all X_dif and X_gcn embeddings
print(f"\n=== Extracting embeddings ===")
embedding_keys = []
for key in adata.obsm.keys():
    if key.startswith('X_dif_nsteps') or key.startswith('X_gcn_nlayers'):
        embedding_keys.append(key)
        print(f"  Found {key}: shape {adata.obsm[key].shape}")

if len(embedding_keys) == 0:
    print(f"Error: No embeddings found matching patterns 'X_dif_nsteps*' or 'X_gcn_nlayers*'")
    print(f"Available obsm keys: {list(adata.obsm.keys())}")
    print("This may indicate that aggregate_embeddings.py did not correctly extract embeddings.")
    sys.exit(1)

# Sort embedding keys for consistent ordering
embedding_keys = sorted(embedding_keys)

# Determine pre-integrated embedding key (use X_fae if available)
pre_integrated_key = None
try:
    if 'X_fae' in adata.obsm:
        pre_integrated_key = 'X_fae'
        print(f"\nUsing 'X_fae' as pre-integrated embedding")
    elif 'X_pca' in adata.obsm:
        pre_integrated_key = 'X_pca'
        print(f"\nUsing 'X_pca' as pre-integrated embedding")
    else:
        print("\nWarning: No suitable pre-integrated embedding found, will skip PCR comparison")
except Exception as e:
    print(f"Warning: Error checking for pre-integrated embedding: {e}")

# Set up batch correction metrics
try:
    batch_corr = BatchCorrection(pcr_comparison=(pre_integrated_key is not None))
    print("BatchCorrection metrics configured successfully")
except Exception as e:
    print(f"Error creating BatchCorrection: {e}")
    traceback.print_exc()
    sys.exit(1)

print(f"\n=== Running SCIB benchmark ===")
print(f"  Batch key: {params.batch_key}")
print(f"  Label key: {params.label_key}")
print(f"  Number of embeddings: {len(embedding_keys)}")
if pre_integrated_key:
    print(f"  Pre-integrated key: {pre_integrated_key}")

# Create Benchmarker
try:
    print(f"\n=== Creating Benchmarker ===")
    bm = Benchmarker(
        adata,
        batch_key=params.batch_key,
        label_key=params.label_key,
        embedding_obsm_keys=embedding_keys,
        pre_integrated_embedding_obsm_key=pre_integrated_key,
        batch_correction_metrics=batch_corr,
        n_jobs=params.n_jobs if hasattr(params, 'n_jobs') else 1,
    )
    print("Benchmarker created successfully")
except Exception as e:
    print(f"Error creating Benchmarker: {e}")
    print(f"  batch_key: {params.batch_key}")
    print(f"  label_key: {params.label_key}")
    print(f"  embedding_keys: {embedding_keys}")
    print(f"  pre_integrated_key: {pre_integrated_key}")
    traceback.print_exc()
    sys.exit(1)

# Run benchmark
try:
    print("\nStarting benchmark...")
    bm.benchmark()
    print("Benchmark completed successfully")
except Exception as e:
    print(f"Error running benchmark: {e}")
    traceback.print_exc()
    sys.exit(1)

print("\n=== Benchmark complete ===")

# Save results
try:
    print(f"\n=== Saving results ===")

    # Save results table as CSV
    results_df = bm.get_results(min_max_scale=False)
    print(f"Saving results table to: {output_table}")
    Path(output_table).parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_table, index=True)
    print(f"Results table shape: {results_df.shape}")
    print("\nResults summary:")
    print(results_df)
except Exception as e:
    print(f"Error saving results table: {e}")
    traceback.print_exc()
    sys.exit(1)

# Generate and save plots
print(f"\n=== Generating and saving plots ===")
print(f"Generating results table plot and saving to: {output_plot}")
Path(output_plot).parent.mkdir(parents=True, exist_ok=True)

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

# Ensure output file exists (create placeholder if plotting failed)
if not Path(output_plot).exists():
    print(f"  Creating placeholder plot file...")
    Path(output_plot).parent.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, 'Plot generation failed.\nCheck logs for details.', 
                ha='center', va='center', fontsize=12, 
                transform=ax.transAxes)
        ax.set_title('SCIB Results Plot - Generation Failed')
        plt.savefig(output_plot, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Created placeholder plot file: {output_plot}")
    except Exception as e:
        Path(output_plot).touch()
        print(f"  Warning: Could not create placeholder plot: {e}")

# Save full benchmarker object (pickle)
try:
    print(f"\nSaving benchmarker object to: {output_results}")
    Path(output_results).parent.mkdir(parents=True, exist_ok=True)
    with open(output_results, 'wb') as f:
        pickle.dump(bm, f)
    print("Benchmarker object saved successfully")
except Exception as e:
    print(f"Error saving benchmarker object: {e}")
    traceback.print_exc()
    sys.exit(1)

print("\n=== SCIB evaluation complete! ===")

