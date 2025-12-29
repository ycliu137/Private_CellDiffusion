"""
Plot combined UMAP visualizations: CellDiffusion and GCN for each network layer
Each page shows both methods (batch and labels) for the same network layer
"""
import sys
from pathlib import Path
import re

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Get input files and output from snakemake
try:
    input_h5ad_files = snakemake.input.h5ad_files
    output_pdf = snakemake.output.pdf
    params = snakemake.params
except Exception as e:
    print(f"Error getting snakemake inputs/outputs: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print(f"\n=== Loading UMAP data from {len(input_h5ad_files)} files ===")

# Organize files by network layer
layer_data = {}  # {layer_number: {'celldiffusion': adata_path, 'gcn': adata_path}}

for h5ad_file in input_h5ad_files:
    filename = Path(h5ad_file).name
    
    # Pattern: celldiffusion_umap_nsteps{N}.h5ad or gcn_umap_nlayers{N}.h5ad
    if 'celldiffusion_umap_nsteps' in filename:
        match = re.search(r'celldiffusion_umap_nsteps(\d+)\.h5ad', filename)
        if match:
            layer = int(match.group(1))
            if layer not in layer_data:
                layer_data[layer] = {}
            layer_data[layer]['celldiffusion'] = h5ad_file
            print(f"  Found CellDiffusion layer {layer}: {filename}")
    elif 'gcn_umap_nlayers' in filename:
        match = re.search(r'gcn_umap_nlayers(\d+)\.h5ad', filename)
        if match:
            layer = int(match.group(1))
            if layer not in layer_data:
                layer_data[layer] = {}
            layer_data[layer]['gcn'] = h5ad_file
            print(f"  Found GCN layer {layer}: {filename}")

# Sort layers
sorted_layers = sorted(layer_data.keys())
print(f"\n=== Found {len(sorted_layers)} network layers: {sorted_layers} ===")

# Check if we have any data
if len(sorted_layers) == 0:
    print("Error: No UMAP data files found!")
    print(f"Expected files matching patterns:")
    print(f"  - celldiffusion_umap_nsteps*.h5ad")
    print(f"  - gcn_umap_nlayers*.h5ad")
    sys.exit(1)

# Load data for each layer (we'll load on demand for each page)
print(f"\n=== Creating multi-page PDF ===")
print(f"Output file: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)

# Create PDF with multiple pages
with PdfPages(output_pdf) as pdf:
    for layer in sorted_layers:
        print(f"\n  Processing layer {layer}...")
        
        # Create figure with 2x2 subplots: 
        # Row 1: CellDiffusion (batch, labels)
        # Row 2: GCN (batch, labels)
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))
        fig.suptitle(f'Network Layers: {layer}', fontsize=16, fontweight='bold', y=0.98)
        
        # Load CellDiffusion data if available
        if 'celldiffusion' in layer_data[layer]:
            try:
                print(f"    Loading CellDiffusion data for layer {layer}...")
                adata_dif = sc.read_h5ad(layer_data[layer]['celldiffusion'])
            umap_key_dif = f'X_umap_dif_nsteps{layer}'
            
            if umap_key_dif in adata_dif.obsm:
                # Temporarily set UMAP for plotting
                adata_dif.obsm['X_umap'] = adata_dif.obsm[umap_key_dif].copy()
                
                # Plot CellDiffusion - Batch (top left)
                sc.pl.umap(
                    adata_dif,
                    color=params.batch_key,
                    ax=axes[0, 0],
                    show=False,
                    frameon=False,
                    title=f"CellDiffusion - Batch"
                )
                
                # Plot CellDiffusion - Labels (top right)
                sc.pl.umap(
                    adata_dif,
                    color=params.label_key,
                    ax=axes[0, 1],
                    show=False,
                    frameon=False,
                    title=f"CellDiffusion - Labels"
                )
                
                # Clean up
                del adata_dif.obsm['X_umap']
            else:
                print(f"    Warning: {umap_key_dif} not found in CellDiffusion data")
                axes[0, 0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
                axes[0, 1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            except Exception as e:
                print(f"    Error loading CellDiffusion data: {e}")
                axes[0, 0].text(0.5, 0.5, 'Error loading data', ha='center', va='center')
                axes[0, 1].text(0.5, 0.5, 'Error loading data', ha='center', va='center')
        else:
            print(f"    Warning: CellDiffusion data not found for layer {layer}")
            axes[0, 0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            axes[0, 1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        
        # Load GCN data if available
        if 'gcn' in layer_data[layer]:
            try:
                print(f"    Loading GCN data for layer {layer}...")
                adata_gcn = sc.read_h5ad(layer_data[layer]['gcn'])
            umap_key_gcn = f'X_umap_gcn_nlayers{layer}'
            
            if umap_key_gcn in adata_gcn.obsm:
                # Temporarily set UMAP for plotting
                adata_gcn.obsm['X_umap'] = adata_gcn.obsm[umap_key_gcn].copy()
                
                # Plot GCN - Batch (bottom left)
                sc.pl.umap(
                    adata_gcn,
                    color=params.batch_key,
                    ax=axes[1, 0],
                    show=False,
                    frameon=False,
                    title=f"GCN - Batch"
                )
                
                # Plot GCN - Labels (bottom right)
                sc.pl.umap(
                    adata_gcn,
                    color=params.label_key,
                    ax=axes[1, 1],
                    show=False,
                    frameon=False,
                    title=f"GCN - Labels"
                )
                
                # Clean up
                del adata_gcn.obsm['X_umap']
            else:
                print(f"    Warning: {umap_key_gcn} not found in GCN data")
                axes[1, 0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
                axes[1, 1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            except Exception as e:
                print(f"    Error loading GCN data: {e}")
                axes[1, 0].text(0.5, 0.5, 'Error loading data', ha='center', va='center')
                axes[1, 1].text(0.5, 0.5, 'Error loading data', ha='center', va='center')
        else:
            print(f"    Warning: GCN data not found for layer {layer}")
            axes[1, 0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            axes[1, 1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for suptitle
        pdf.savefig(fig, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved page for layer {layer}")

print(f"\n=== Multi-page PDF created successfully! ===")
print(f"  Total pages: {len(sorted_layers)}")
print(f"  Output file: {output_pdf}")

