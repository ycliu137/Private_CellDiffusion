"""
SCIB evaluation for BBKNN integration (single embedding X_bbknn).
"""
import sys
from pathlib import Path
import pickle
import traceback

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scanpy as sc

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))

from scib_metrics.benchmark import BatchCorrection, Benchmarker

input_h5ad = snakemake.input.h5ad
output_results = snakemake.output.results
output_table = snakemake.output.table
output_plot = snakemake.output.plot
params = snakemake.params

print("\n=== Loading integrated adata ===")
print(f"Input: {input_h5ad}")
adata = sc.read_h5ad(input_h5ad)
print(f"Shape: {adata.shape}")

embedding_key = "X_bbknn"
if embedding_key not in adata.obsm:
    raise ValueError(f"obsm['{embedding_key}'] not found. Keys: {list(adata.obsm.keys())}")
all_embedding_keys = [embedding_key]

pre_integrated_key = None
if "X_pca" in adata.obsm:
    pre_integrated_key = "X_pca"
    print("Using 'X_pca' as pre-integrated embedding")
else:
    cand = [k for k in adata.obsm.keys() if k not in (embedding_key, "latent")]
    if cand:
        pre_integrated_key = cand[0]
        print(f"Using '{pre_integrated_key}' as pre-integrated (fallback)")
    else:
        print("No pre-integrated embedding; PCR comparison skipped")

batch_corr = BatchCorrection(pcr_comparison=False)
print("\n=== Running SCIB benchmark ===")
print(f"  Batch key: {params.batch_key}")
print(f"  Label key: {params.label_key}")
print(f"  Embedding: {all_embedding_keys}")

bm = Benchmarker(
    adata,
    batch_key=params.batch_key,
    label_key=params.label_key,
    embedding_obsm_keys=all_embedding_keys,
    pre_integrated_embedding_obsm_key=pre_integrated_key,
    batch_correction_metrics=batch_corr,
    n_jobs=getattr(params, "n_jobs", 1),
)
bm.benchmark()
print("Benchmark complete.")

Path(output_table).parent.mkdir(parents=True, exist_ok=True)
results_df = bm.get_results(min_max_scale=False)
results_df.to_csv(output_table, index=True)
print(f"Results table saved: {output_table}")

plot_saved = False
try:
    try:
        from scib_metrics_new.benchmarked_plot_results_table import plot_results_table_new
        bm.plot_results_table_new = plot_results_table_new
        fig = bm.plot_results_table_new(min_max_scale=False, show=False, save_fig=str(output_plot))
        if fig is not None:
            import matplotlib.pyplot as plt
            plt.close(fig)
        if Path(output_plot).exists():
            plot_saved = True
    except Exception:
        pass
    if not plot_saved:
        import matplotlib.pyplot as plt
        fig = bm.plot_results_table(min_max_scale=False)
        if fig is not None:
            plt.savefig(output_plot, dpi=300, bbox_inches="tight")
            plt.close(fig)
            plot_saved = Path(output_plot).exists()
except Exception as e:
    print(f"Plot failed: {e}")
    traceback.print_exc()

if not Path(output_plot).exists():
    Path(output_plot).parent.mkdir(parents=True, exist_ok=True)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, "Plot generation failed.\nCheck logs.", ha="center", va="center", transform=ax.transAxes)
    ax.axis("off")
    plt.savefig(output_plot, dpi=150, bbox_inches="tight")
    plt.close(fig)

Path(output_results).parent.mkdir(parents=True, exist_ok=True)
with open(output_results, "wb") as f:
    pickle.dump(bm, f)
print(f"Benchmarker saved: {output_results}")
print("=== SCIB evaluation complete ===")
