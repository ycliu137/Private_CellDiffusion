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

# Create figure with subplots: 2 rows (Batch, Labels) x 2 cols (CellDiffusion, UniPort)
fig = plt.figure(figsize=(14, 12))
gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3,
              left=0.08, right=0.95, top=0.95, bottom=0.1)

axes = [[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(2)]

# Set equal aspect ratio for square plots
for i in range(2):
    for j in range(2):
        axes[i][j].set_aspect('equal', adjustable='box')

print(f"\n=== Plotting UMAPs ===")

# Row 0: Batch coloring
# Column 0: CellDiffusion
print(f"  Plotting CellDiffusion Batch UMAP...")
if cdif_umap_key != 'X_umap':
    adata_cdif.obsm['X_umap'] = adata_cdif.obsm[cdif_umap_key].copy()
sc.pl.umap(
    adata_cdif,
    color=batch_key,
    ax=axes[0][0],
    show=False,
    frameon=False,
    title="CellDiffusion - Batch",
    legend_loc='none'
)
axes[0][0].title.set_fontsize(14)
if cdif_umap_key != 'X_umap':
    del adata_cdif.obsm['X_umap']

# Column 1: UniPort
print(f"  Plotting UniPort Batch UMAP...")
if uniport_umap_key != 'X_umap':
    adata_uniport.obsm['X_umap'] = adata_uniport.obsm[uniport_umap_key].copy()
sc.pl.umap(
    adata_uniport,
    color=batch_key,
    ax=axes[0][1],
    show=False,
    frameon=False,
    title="UniPort - Batch",
    legend_loc='none'
)
axes[0][1].title.set_fontsize(14)
if uniport_umap_key != 'X_umap':
    del adata_uniport.obsm['X_umap']

# Row 1: Cell Type coloring
# Column 0: CellDiffusion
print(f"  Plotting CellDiffusion Labels UMAP...")
if cdif_umap_key != 'X_umap':
    adata_cdif.obsm['X_umap'] = adata_cdif.obsm[cdif_umap_key].copy()
sc.pl.umap(
    adata_cdif,
    color=label_key,
    ax=axes[1][0],
    show=False,
    frameon=False,
    title="CellDiffusion - Cell Type",
    legend_loc='none'
)
axes[1][0].title.set_fontsize(14)
if cdif_umap_key != 'X_umap':
    del adata_cdif.obsm['X_umap']

# Column 1: UniPort
print(f"  Plotting UniPort Labels UMAP...")
if uniport_umap_key != 'X_umap':
    adata_uniport.obsm['X_umap'] = adata_uniport.obsm[uniport_umap_key].copy()
sc.pl.umap(
    adata_uniport,
    color=label_key,
    ax=axes[1][1],
    show=False,
    frameon=False,
    title="UniPort - Cell Type",
    legend_loc='none'
)
axes[1][1].title.set_fontsize(14)
if uniport_umap_key != 'X_umap':
    del adata_uniport.obsm['X_umap']

# Extract and add legends
print(f"\n=== Extracting legends ===")

# Extract batch legend from CellDiffusion
if cdif_umap_key != 'X_umap':
    adata_cdif.obsm['X_umap'] = adata_cdif.obsm[cdif_umap_key].copy()
temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
sc.pl.umap(adata_cdif, color=batch_key, ax=temp_ax, show=False, frameon=False)
batch_handles, batch_labels = temp_ax.get_legend_handles_labels()
plt.close(temp_fig)
if cdif_umap_key != 'X_umap':
    del adata_cdif.obsm['X_umap']

# Extract label legend from CellDiffusion
if cdif_umap_key != 'X_umap':
    adata_cdif.obsm['X_umap'] = adata_cdif.obsm[cdif_umap_key].copy()
temp_fig, temp_ax = plt.subplots(figsize=(1, 1))
sc.pl.umap(adata_cdif, color=label_key, ax=temp_ax, show=False, frameon=False)
label_handles, label_labels = temp_ax.get_legend_handles_labels()
plt.close(temp_fig)
if cdif_umap_key != 'X_umap':
    del adata_cdif.obsm['X_umap']

# Add legends below the plots
fig.text(0.25, 0.02, 'Batch Legend:', fontsize=12, fontweight='bold', ha='center')
fig.text(0.75, 0.02, 'Cell Type Legend:', fontsize=12, fontweight='bold', ha='center')

# Create legend items (simplified for now - just use text)
if batch_handles and batch_labels:
    batch_text = ', '.join(batch_labels[:5])  # Show first 5
    if len(batch_labels) > 5:
        batch_text += f', ... ({len(batch_labels)} total)'
    # fig.text(0.25, -0.02, batch_text, fontsize=10, ha='center')

if label_handles and label_labels:
    label_text = ', '.join(label_labels[:5])  # Show first 5
    if len(label_labels) > 5:
        label_text += f', ... ({len(label_labels)} total)'
    # fig.text(0.75, -0.02, label_text, fontsize=10, ha='center')

# Save figure
print(f"\n=== Saving comparison plot ===")
print(f"  Output: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print("\n=== Comparison plot complete! ===")
