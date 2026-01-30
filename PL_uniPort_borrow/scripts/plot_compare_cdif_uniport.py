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
import numpy as np

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

# Ensure categorical variables and sync their color mapping
print(f"\n=== Setting up unified color mapping ===")

# Convert to categorical if not already
for adata in [adata_cdif, adata_uniport]:
    if batch_key not in adata.obs.columns:
        print(f"ERROR: {batch_key} not found in data")
        sys.exit(1)
    if not isinstance(adata.obs[batch_key], type(adata_cdif.obs[batch_key])):
        if not adata.obs[batch_key].dtype.name == 'category':
            adata.obs[batch_key] = adata.obs[batch_key].astype('category')
    
    if hasattr(params, 'label_key') and params.label_key and label_key in adata.obs.columns:
        if not adata.obs[label_key].dtype.name == 'category':
            adata.obs[label_key] = adata.obs[label_key].astype('category')
    
    if has_lineage and 'lineage' in adata.obs.columns:
        if not adata.obs['lineage'].dtype.name == 'category':
            adata.obs['lineage'] = adata.obs['lineage'].astype('category')

# Create unified categories from CellDiffusion data
batch_categories = list(adata_cdif.obs[batch_key].cat.categories)
if label_key in adata_cdif.obs.columns:
    label_categories = list(adata_cdif.obs[label_key].cat.categories)
else:
    label_categories = None

if has_lineage:
    lineage_categories = list(adata_cdif.obs['lineage'].cat.categories)
else:
    lineage_categories = None

# Sync categories in both datasets
for adata in [adata_cdif, adata_uniport]:
    adata.obs[batch_key] = adata.obs[batch_key].cat.set_categories(batch_categories)
    
    if label_categories and label_key in adata.obs.columns:
        adata.obs[label_key] = adata.obs[label_key].cat.set_categories(label_categories)
    
    if lineage_categories and 'lineage' in adata.obs.columns:
        adata.obs['lineage'] = adata.obs['lineage'].cat.set_categories(lineage_categories)

# Generate color mapping using scanpy's color schema and store in .uns
print(f"  Generating color mapping...")
from matplotlib import cm
from matplotlib.colors import rgb2hex

# Generate palette for batch
batch_palette = sc.pl.palettes.default_20
batch_colors = [rgb2hex(batch_palette[i % len(batch_palette)]) for i in range(len(batch_categories))]
adata_cdif.uns[f'{batch_key}_colors'] = batch_colors
adata_uniport.uns[f'{batch_key}_colors'] = batch_colors

print(f"  Batch categories: {batch_categories}")
print(f"  Batch colors assigned: {len(batch_colors)}")

# Generate palette for labels
if label_categories:
    label_colors = [rgb2hex(batch_palette[i % len(batch_palette)]) for i in range(len(label_categories))]
    adata_cdif.uns[f'{label_key}_colors'] = label_colors
    adata_uniport.uns[f'{label_key}_colors'] = label_colors
    print(f"  Label categories: {label_categories}")
    print(f"  Label colors assigned: {len(label_colors)}")

# Generate palette for lineage
if lineage_categories:
    lineage_colors = [rgb2hex(batch_palette[i % len(batch_palette)]) for i in range(len(lineage_categories))]
    adata_cdif.uns[f'lineage_colors'] = lineage_colors
    adata_uniport.uns[f'lineage_colors'] = lineage_colors
    print(f"  Lineage categories: {lineage_categories}")
    print(f"  Lineage colors assigned: {len(lineage_colors)}")

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
    if hasattr(params, 'label_key') and params.label_key and label_key in madata.obs.columns:
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
    if has_lineage and 'lineage' in madata.obs.columns:
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

# Extract legends from the color mapping in .uns
print(f"\n=== Extracting legends ===")
if cdif_umap_key != 'X_umap':
    adata_cdif.obsm['X_umap'] = adata_cdif.obsm[cdif_umap_key].copy()

temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
sc.pl.umap(adata_cdif, color=batch_key, ax=temp_ax, show=False, frameon=False)
batch_handles, batch_labels = temp_ax.get_legend_handles_labels()
plt.close(temp_fig)

if label_categories:
    temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
    sc.pl.umap(adata_cdif, color=label_key, ax=temp_ax, show=False, frameon=False)
    label_handles, label_labels = temp_ax.get_legend_handles_labels()
    plt.close(temp_fig)
else:
    label_handles, label_labels = None, None

if lineage_categories:
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
