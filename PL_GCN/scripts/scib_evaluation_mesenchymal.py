"""
SCIB evaluation for Mesenchymal lineage cells only: benchmark CellDiffusion (X_dif) vs GCN (X_gcn) embeddings
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
from scib_metrics.benchmark import BioConservation

# Get input and output from snakemake
input_h5ad = snakemake.input.h5ad
output_results = snakemake.output.results
output_table = snakemake.output.table
output_plot = snakemake.output.plot
params = snakemake.params

print(f"\n=== Loading adata ===")
print(f"Input file: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)
print(f"Original data shape: {adata.shape}")

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

# Get X_dif and X_gcn embeddings
print(f"\n=== Extracting embeddings ===")
embedding_keys = []
if 'X_dif' in adata.obsm:
    embedding_keys.append('X_dif')
    print(f"  Found X_dif: shape {adata.obsm['X_dif'].shape}")
else:
    print(f"  Warning: X_dif not found in adata.obsm")

if 'X_gcn' in adata.obsm:
    embedding_keys.append('X_gcn')
    print(f"  Found X_gcn: shape {adata.obsm['X_gcn'].shape}")
else:
    print(f"  Warning: X_gcn not found in adata.obsm")

if len(embedding_keys) == 0:
    raise ValueError("No embeddings found (X_dif or X_gcn) in adata.obsm")

# Determine pre-integrated embedding key (use X_fae if available)
pre_integrated_key = None
if 'X_fae' in adata.obsm:
    pre_integrated_key = 'X_fae'
    print(f"\nUsing 'X_fae' as pre-integrated embedding")
elif 'X_pca' in adata.obsm:
    pre_integrated_key = 'X_pca'
    print(f"\nUsing 'X_pca' as pre-integrated embedding")
else:
    print("\nWarning: No suitable pre-integrated embedding found")

# Set up batch correction metrics (PCR comparison disabled)
batch_corr = BatchCorrection(pcr_comparison=False)
bio_cons = BioConservation(isolated_labels=False, clisi_knn=False)

print(f"\n=== Running SCIB benchmark (Mesenchymal lineage only) ===")
print(f"  Batch key: {params.batch_key}")
print(f"  Label key: {params.label_key}")
print(f"  Embedding keys: {embedding_keys}")
if pre_integrated_key:
    print(f"  Pre-integrated key: {pre_integrated_key}")

# Create Benchmarker
bm = Benchmarker(
    adata,
    batch_key=params.batch_key,
    label_key=params.label_key,
    embedding_obsm_keys=embedding_keys,
    pre_integrated_embedding_obsm_key=pre_integrated_key,
    batch_correction_metrics=batch_corr,
    bio_conservation_metrics=bio_cons,
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
        
        # Verify file was created
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
                # Verify file was created
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
# This prevents Snakemake from failing due to missing output file
if not Path(output_plot).exists():
    print(f"  Creating empty plot file as placeholder...")
    Path(output_plot).parent.mkdir(parents=True, exist_ok=True)
    # Create a simple placeholder PDF file
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, 'Plot generation failed.\nCheck logs for details.', 
                ha='center', va='center', fontsize=12, 
                transform=ax.transAxes)
        ax.set_title('SCIB Results Plot (Mesenchymal) - Generation Failed')
        plt.savefig(output_plot, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Created placeholder plot file: {output_plot}")
    except Exception as e:
        # If even placeholder fails, create empty file
        Path(output_plot).touch()
        print(f"  Warning: Could not create placeholder plot: {e}")
        print(f"  Created empty file as last resort: {output_plot}")

# Save full benchmarker object (pickle)
print(f"\nSaving benchmarker object to: {output_results}")
Path(output_results).parent.mkdir(parents=True, exist_ok=True)
with open(output_results, 'wb') as f:
    pickle.dump(bm, f)

print("\n=== SCIB evaluation (Mesenchymal) complete! ===")

