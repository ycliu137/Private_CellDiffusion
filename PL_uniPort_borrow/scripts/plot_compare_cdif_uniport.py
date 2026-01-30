"""
Plot UMAP comparison between CellDiffusion and UniPort integration
Displays side-by-side: Batch and Cell Type UMAPs
"""
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec

# Set up figure parameters
rcParams["figure.figsize"] = (12, 6)
rcParams["font.size"] = 12
rcParams["axes.titlesize"] = 14
rcParams["axes.labelsize"] = 12

# Get inputs from snakemake
cdif_h5ad = snakemake.input.cdif_h5ad
uniport_h5ad = snakemake.input.uniport_h5ad
output_pdf = snakemake.output.pdf
params = snakemake.params

batch_key = params.batch_key
label_key = params.label_key

print(f"\n=== Loading CellDiffusion data ===")
print(f"  File: {cdif_h5ad}")
adata_cdif = sc.read_h5ad(cdif_h5ad)
print(f"  Shape: {adata_cdif.shape}")
print(f"  Available obsm keys: {list(adata_cdif.obsm.keys())}")

print(f"\n=== Loading UniPort data ===")
print(f"  File: {uniport_h5ad}")
adata_uniport = sc.read_h5ad(uniport_h5ad)
print(f"  Shape: {adata_uniport.shape}")
print(f"  Available obsm keys: {list(adata_uniport.obsm.keys())}")

# Check for UMAP embeddings
cdif_umap_key = None
uniport_umap_key = None

# Find CellDiffusion UMAP key (X_umap_dif)
if 'X_umap_dif' in adata_cdif.obsm.keys():
    cdif_umap_key = 'X_umap_dif'
    print(f"  Found CellDiffusion UMAP: {cdif_umap_key}")
elif 'X_umap' in adata_cdif.obsm.keys():
    cdif_umap_key = 'X_umap'
    print(f"  Found UMAP: {cdif_umap_key}")
else:
    # Find any umap key
    for key in adata_cdif.obsm.keys():
        if 'umap' in key.lower():
            cdif_umap_key = key
            print(f"  Found UMAP: {cdif_umap_key}")
            break

# Find UniPort UMAP key
if 'X_umap' in adata_uniport.obsm.keys():
    uniport_umap_key = 'X_umap'
    print(f"  Found UniPort UMAP: {uniport_umap_key}")
else:
    # Find any umap key
    for key in adata_uniport.obsm.keys():
        if 'umap' in key.lower():
            uniport_umap_key = key
            print(f"  Found UniPort UMAP: {uniport_umap_key}")
            break

if not cdif_umap_key:
    print(f"ERROR: No UMAP found in CellDiffusion data")
    sys.exit(1)

if not uniport_umap_key:
    print(f"ERROR: No UMAP found in UniPort data. Computing UMAP...")
    # Compute UMAP if not present
    if 'X_uniport' in adata_uniport.obsm.keys():
        sc.pp.neighbors(adata_uniport, use_rep='X_uniport')
        sc.tl.umap(adata_uniport)
        uniport_umap_key = 'X_umap'
    else:
        print(f"ERROR: Cannot compute UMAP - no embedding found")
        sys.exit(1)

# Create figure with subplots: 2 methods (CellDiffusion, UniPort) x 3 cols (Batch, Labels, Lineage)
n_methods = 2
n_cols = 3
subplot_size = 5
fig_width = subplot_size * n_cols + 2
fig_height = subplot_size * n_methods + 1.5

fig = plt.figure(figsize=(fig_width, fig_height))
gs = GridSpec(n_methods + 1, n_cols, figure=fig,
              hspace=0.3, wspace=0.3,
              left=0.08, right=0.95, top=0.95, bottom=0.08)

# Create axes array (last row reserved for legends)
axes = []
for i in range(n_methods + 1):
    row_axes = []
    for j in range(n_cols):
        ax = fig.add_subplot(gs[i, j])
        if i < n_methods:
            ax.set_aspect('equal', adjustable='box')
        row_axes.append(ax)
    axes.append(row_axes)

has_lineage = 'lineage' in adata_cdif.obs.columns

print(f"\n=== Plotting UMAPs ===")

methods = [("CellDiffusion", adata_cdif, cdif_umap_key), ("UniPort", adata_uniport, uniport_umap_key)]

for i, (mname, madata, mukey) in enumerate(methods):
    # Batch
    print(f"  Plotting {mname} Batch UMAP...")
    if mukey != 'X_umap':
        madata.obsm['X_umap'] = madata.obsm[mukey].copy()
    sc.pl.umap(madata, color=batch_key, ax=axes[i][0], show=False, frameon=False,
               title=f"{mname} - Batch", legend_loc='none')
    axes[i][0].title.set_fontsize(14)
    if mukey != 'X_umap':
        del madata.obsm['X_umap']

    # Labels
    print(f"  Plotting {mname} Labels UMAP...")
    if hasattr(params, 'label_key') and params.label_key:
        if mukey != 'X_umap':
            madata.obsm['X_umap'] = madata.obsm[mukey].copy()
        sc.pl.umap(madata, color=label_key, ax=axes[i][1], show=False, frameon=False,
                   title=f"{mname} - Labels", legend_loc='none')
        axes[i][1].title.set_fontsize(14)
        if mukey != 'X_umap':
            del madata.obsm['X_umap']
    else:
        axes[i][1].axis('off')

    # Lineage
    if has_lineage:
        print(f"  Plotting {mname} Lineage UMAP...")
        if mukey != 'X_umap':
            madata.obsm['X_umap'] = madata.obsm[mukey].copy()
        sc.pl.umap(madata, color='lineage', ax=axes[i][2], show=False, frameon=False,
                   title=f"{mname} - Lineage", legend_loc='none')
        axes[i][2].title.set_fontsize(14)
        if mukey != 'X_umap':
            del madata.obsm['X_umap']
    else:
        axes[i][2].axis('off')

# Extract legends from CellDiffusion (use first method's embedding)
print(f"\n=== Extracting legends ===")
if cdif_umap_key != 'X_umap':
    adata_cdif.obsm['X_umap'] = adata_cdif.obsm[cdif_umap_key].copy()

temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
sc.pl.umap(adata_cdif, color=batch_key, ax=temp_ax, show=False, frameon=False)
batch_handles, batch_labels = temp_ax.get_legend_handles_labels()
plt.close(temp_fig)

if hasattr(params, 'label_key') and params.label_key:
    temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
    sc.pl.umap(adata_cdif, color=label_key, ax=temp_ax, show=False, frameon=False)
    label_handles, label_labels = temp_ax.get_legend_handles_labels()
    plt.close(temp_fig)
else:
    label_handles, label_labels = None, None

if has_lineage:
    temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
    sc.pl.umap(adata_cdif, color='lineage', ax=temp_ax, show=False, frameon=False)
    lineage_handles, lineage_labels = temp_ax.get_legend_handles_labels()
    plt.close(temp_fig)
else:
    lineage_handles, lineage_labels = None, None

if cdif_umap_key != 'X_umap':
    del adata_cdif.obsm['X_umap']

# Place legends in the bottom (legend) row
legend_y_top = 0.92

# Column 0: Batch legend
if batch_handles and batch_labels:
    axes[n_methods][0].axis('off')
    axes[n_methods][0].text(0.5, 0.98, 'Batch', transform=axes[n_methods][0].transAxes,
                            ha='center', va='top', fontsize=12, fontweight='bold')
    axes[n_methods][0].legend(batch_handles, batch_labels, loc='upper center',
                              frameon=False, fontsize=10, ncol=1, bbox_to_anchor=(0.5, legend_y_top))
else:
    axes[n_methods][0].axis('off')

# Column 1: Labels legend
if label_handles and label_labels:
    axes[n_methods][1].axis('off')
    axes[n_methods][1].text(0.5, 0.98, 'Labels', transform=axes[n_methods][1].transAxes,
                            ha='center', va='top', fontsize=12, fontweight='bold')
    axes[n_methods][1].legend(label_handles, label_labels, loc='upper center', frameon=False,
                              fontsize=10, ncol=3, bbox_to_anchor=(0.5, legend_y_top))
else:
    axes[n_methods][1].axis('off')

# Column 2: Lineage legend
if lineage_handles and lineage_labels:
    axes[n_methods][2].axis('off')
    axes[n_methods][2].text(0.5, 0.98, 'Lineage', transform=axes[n_methods][2].transAxes,
                            ha='center', va='top', fontsize=12, fontweight='bold')
    axes[n_methods][2].legend(lineage_handles, lineage_labels, loc='upper center', frameon=False,
                              fontsize=10, ncol=1, bbox_to_anchor=(0.5, legend_y_top))
else:
    axes[n_methods][2].axis('off')

# Save figure
print(f"\n=== Saving comparison plot ===")
print(f"  Output: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print("\n=== Comparison plot complete! ===")
