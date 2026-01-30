"""
Compare CellDiffusion and BBKNN integration results
"""
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Load outputs
celldiffusion_h5ad = snakemake.input.celldiffusion_h5ad
bbknn_h5ad = snakemake.input.bbknn_h5ad
output_pdf = snakemake.output.pdf
params = snakemake.params

print(f"\n=== Loading integrated adata ===")
print(f"  CellDiffusion: {celldiffusion_h5ad}")
print(f"  BBKNN: {bbknn_h5ad}")

adata_cdif = sc.read_h5ad(celldiffusion_h5ad)
adata_bbknn = sc.read_h5ad(bbknn_h5ad)

print(f"  CellDiffusion shape: {adata_cdif.shape}")
print(f"  BBKNN shape: {adata_bbknn.shape}")

# Ensure same cell order for color mapping
common_cells = list(set(adata_cdif.obs_names) & set(adata_bbknn.obs_names))
adata_cdif = adata_cdif[common_cells]
adata_bbknn = adata_bbknn[common_cells]

print(f"  Common cells: {len(common_cells)}")

# Create comparison plot: 1x4 grid
# [CellDif Batch] [CellDif Labels] [BBKNN Batch] [BBKNN Labels]
fig = plt.figure(figsize=(16, 4))
gs = fig.add_gridspec(2, 4, height_ratios=[20, 1], hspace=0.3, wspace=0.3)

# CellDiffusion - Batch
ax1 = fig.add_subplot(gs[0, 0])
sc.pl.umap(adata_cdif, color=params.batch_key, ax=ax1, show=False, 
           title='CellDiffusion - Batch', legend_loc=None, frameon=False)

# CellDiffusion - Labels
ax2 = fig.add_subplot(gs[0, 1])
sc.pl.umap(adata_cdif, color=params.label_key, ax=ax2, show=False,
           title='CellDiffusion - Labels', legend_loc=None, frameon=False)

# BBKNN - Batch
ax3 = fig.add_subplot(gs[0, 2])
sc.pl.umap(adata_bbknn, color=params.batch_key, ax=ax3, show=False,
           title='BBKNN - Batch', legend_loc=None, frameon=False)

# BBKNN - Labels
ax4 = fig.add_subplot(gs[0, 3])
sc.pl.umap(adata_bbknn, color=params.label_key, ax=ax4, show=False,
           title='BBKNN - Labels', legend_loc=None, frameon=False)

# Add shared legends in bottom row
# Legend for batch (from CellDiffusion)
batch_colors = adata_cdif.obs[params.batch_key].cat.categories
# Get colors from adata if available, otherwise use default palette
if f"{params.batch_key}_colors" in adata_cdif.uns:
    batch_color_list = adata_cdif.uns[f"{params.batch_key}_colors"]
else:
    import seaborn as sns
    batch_color_list = sns.color_palette("husl", len(batch_colors))
batch_handles = [mpatches.Patch(facecolor=batch_color_list[i], 
                                label=batch_colors[i]) 
                for i in range(len(batch_colors))]

# Legend for labels (from CellDiffusion)
label_colors = adata_cdif.obs[params.label_key].cat.categories
# Get colors from adata if available, otherwise use default palette
if f"{params.label_key}_colors" in adata_cdif.uns:
    label_color_list = adata_cdif.uns[f"{params.label_key}_colors"]
else:
    import seaborn as sns
    label_color_list = sns.color_palette("tab20", len(label_colors))
label_handles = [mpatches.Patch(facecolor=label_color_list[i % len(label_color_list)],
                                label=label_colors[i])
                for i in range(len(label_colors))]

# Add batch legend
ax_batch_leg = fig.add_subplot(gs[1, :2])
ax_batch_leg.legend(handles=batch_handles, loc='center', ncol=5, frameon=False)
ax_batch_leg.set_title('Batch', loc='left', fontsize=10)
ax_batch_leg.axis('off')

# Add label legend
ax_label_leg = fig.add_subplot(gs[1, 2:])
ax_label_leg.legend(handles=label_handles, loc='center', ncol=5, frameon=False)
ax_label_leg.set_title('Cell Type', loc='left', fontsize=10)
ax_label_leg.axis('off')

plt.suptitle('CellDiffusion vs BBKNN Integration Comparison', fontsize=14, y=0.98)

print(f"\n=== Saving comparison plot ===")
print(f"  Output: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=150, bbox_inches='tight')
plt.close()

print("\n=== Comparison plot complete! ===")
