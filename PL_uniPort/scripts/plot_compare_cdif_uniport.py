"""
Plot UMAP comparison between CellDiffusion and UniPort integration
Displays side-by-side: Batch and Cell Type UMAPs (2x2 layout)
"""
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

rcParams["figure.figsize"] = (12, 6)
rcParams["font.size"] = 12
rcParams["axes.titlesize"] = 14
rcParams["axes.labelsize"] = 12

cdif_h5ad = snakemake.input.cdif_h5ad
uniport_h5ad = snakemake.input.uniport_h5ad
output_pdf = snakemake.output.pdf
params = snakemake.params
batch_key = params.batch_key
label_key = params.label_key

print(f"\n=== Loading CellDiffusion data ===")
adata_cdif = sc.read_h5ad(cdif_h5ad)
print(f"  Shape: {adata_cdif.shape}")

print(f"\n=== Loading UniPort data ===")
adata_uniport = sc.read_h5ad(uniport_h5ad)
print(f"  Shape: {adata_uniport.shape}")

cdif_umap_key = "X_umap_dif" if "X_umap_dif" in adata_cdif.obsm else "X_umap"
uniport_umap_key = "X_umap"

print(f"\n=== Aligning cell order ===")
common_cells = set(adata_cdif.obs.index) & set(adata_uniport.obs.index)
common_cells_ordered = [cell for cell in adata_cdif.obs.index if cell in common_cells]
adata_uniport = adata_uniport[common_cells_ordered].copy()
adata_cdif = adata_cdif[common_cells_ordered].copy()
adata_cdif.obsm["X_umap_uniport"] = adata_uniport.obsm[uniport_umap_key].copy()

n_cols = 4
n_rows = 2
fig = plt.figure(figsize=(20, 6))
gs = GridSpec(n_rows, n_cols, figure=fig, hspace=0.25, wspace=0.25,
              left=0.06, right=0.98, top=0.95, bottom=0.12)

axes = [[fig.add_subplot(gs[i, j]) for j in range(n_cols)] for i in range(n_rows)]
for j in range(n_cols):
    axes[0][j].set_aspect("equal", adjustable="box")

if not adata_cdif.obs[batch_key].dtype.name == "category":
    adata_cdif.obs[batch_key] = adata_cdif.obs[batch_key].astype("category")
if label_key in adata_cdif.obs.columns and not adata_cdif.obs[label_key].dtype.name == "category":
    adata_cdif.obs[label_key] = adata_cdif.obs[label_key].astype("category")

print(f"\n=== Plotting UMAPs ===")
methods = [("CellDiffusion", cdif_umap_key), ("UniPort", "X_umap_uniport")]
for idx, (mname, umap_key) in enumerate(methods):
    col_batch = idx * 2
    col_labels = idx * 2 + 1
    adata_cdif.obsm["X_umap"] = adata_cdif.obsm[umap_key].copy()
    sc.pl.umap(adata_cdif, color=batch_key, ax=axes[0][col_batch], show=False, frameon=False,
               title=f"{mname} - Batch", legend_loc="none")
    axes[0][col_batch].title.set_fontsize(14)
    if label_key in adata_cdif.obs.columns:
        sc.pl.umap(adata_cdif, color=label_key, ax=axes[0][col_labels], show=False, frameon=False,
                   title=f"{mname} - Labels", legend_loc="none")
        axes[0][col_labels].title.set_fontsize(14)
    else:
        axes[0][col_labels].axis("off")

del adata_cdif.obsm["X_umap"]

adata_cdif.obsm["X_umap"] = adata_cdif.obsm[cdif_umap_key].copy()
temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
sc.pl.umap(adata_cdif, color=batch_key, ax=temp_ax, show=False)
batch_handles, batch_labels = temp_ax.get_legend_handles_labels()
plt.close(temp_fig)

label_handles, label_labels = None, None
if label_key in adata_cdif.obs.columns:
    temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
    sc.pl.umap(adata_cdif, color=label_key, ax=temp_ax, show=False)
    label_handles, label_labels = temp_ax.get_legend_handles_labels()
    plt.close(temp_fig)

del adata_cdif.obsm["X_umap"]

for j in range(n_cols):
    axes[1][j].axis("off")

if batch_handles:
    axes[1][0].text(0.5, 0.98, "Batch", transform=axes[1][0].transAxes,
                    ha="center", va="top", fontsize=12, fontweight="bold")
    axes[1][0].legend(batch_handles, batch_labels, loc="upper center",
                      frameon=False, fontsize=10, bbox_to_anchor=(0.5, 0.92))

if label_handles:
    axes[1][2].text(0.5, 0.98, "Labels", transform=axes[1][2].transAxes,
                    ha="center", va="top", fontsize=12, fontweight="bold")
    axes[1][2].legend(label_handles, label_labels, loc="upper center", frameon=False,
                      fontsize=10, ncol=3, bbox_to_anchor=(0.5, 0.92))

Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches="tight")
plt.close()

print("\n=== Comparison plot complete! ===")
