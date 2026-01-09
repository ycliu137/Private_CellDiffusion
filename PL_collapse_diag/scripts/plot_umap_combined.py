"""
Plot combined UMAP visualizations: CellDiffusion and GCN for each network layer
Each page shows both methods (batch and labels) for the same network layer
Reads from aggregated_embeddings.h5ad to reduce file I/O
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

# Get input file and output from snakemake
try:
    input_h5ad = snakemake.input.h5ad
    output_pdf = snakemake.output.pdf
    params = snakemake.params
except Exception as e:
    print(f"Error getting snakemake inputs/outputs: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print(f"\n=== Loading aggregated UMAP data ===")
print(f"Input file: {input_h5ad}")

# Load aggregated adata
adata = sc.read_h5ad(input_h5ad)
print(f"Data shape: {adata.shape}")

# Find all UMAP keys
print(f"\n=== Finding UMAP embeddings ===")
umap_keys = [k for k in adata.obsm.keys() if 'umap' in k.lower()]
print(f"Found {len(umap_keys)} UMAP embeddings:")
for key in umap_keys:
    print(f"  - {key}")

if len(umap_keys) == 0:
    raise ValueError("No UMAP embeddings found in aggregated_embeddings.h5ad")

# Organize UMAP keys by network layer
layer_data = {}  # {layer_number: {'celldiffusion': umap_key, 'gcn': umap_key}}

for umap_key in umap_keys:
    # Pattern: X_umap_dif_nsteps{N} or X_umap_gcn_nlayers{N}
    if 'X_umap_dif_nsteps' in umap_key:
        match = re.search(r'X_umap_dif_nsteps(\d+)', umap_key)
        if match:
            layer = int(match.group(1))
            if layer not in layer_data:
                layer_data[layer] = {}
            layer_data[layer]['celldiffusion'] = umap_key
            print(f"  Found CellDiffusion layer {layer}: {umap_key}")
    elif 'X_umap_gcn_nlayers' in umap_key:
        match = re.search(r'X_umap_gcn_nlayers(\d+)', umap_key)
        if match:
            layer = int(match.group(1))
            if layer not in layer_data:
                layer_data[layer] = {}
            layer_data[layer]['gcn'] = umap_key
            print(f"  Found GCN layer {layer}: {umap_key}")

# Sort layers
sorted_layers = sorted(layer_data.keys())
print(f"\n=== Found {len(sorted_layers)} network layers: {sorted_layers} ===")

# Check if we have any data
if len(sorted_layers) == 0:
    raise ValueError("No UMAP data found matching expected patterns (X_umap_dif_nsteps* or X_umap_gcn_nlayers*)")

# Create PDF with multiple pages
print(f"\n=== Creating multi-page PDF ===")
print(f"Output file: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)

with PdfPages(output_pdf) as pdf:
    for layer in sorted_layers:
        print(f"\n  Processing layer {layer}...")
        
        # Create figure with 2x2 subplots: 
        # Row 1: CellDiffusion (batch, labels)
        # Row 2: GCN (batch, labels)
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))
        fig.suptitle(f'Network Layers: {layer}', fontsize=16, fontweight='bold', y=0.98)
        
        # Plot CellDiffusion if available
        if 'celldiffusion' in layer_data[layer]:
            try:
                umap_key_dif = layer_data[layer]['celldiffusion']
                print(f"    Plotting CellDiffusion using {umap_key_dif}...")
                
                if umap_key_dif in adata.obsm:
                    # Temporarily set UMAP for plotting
                    adata.obsm['X_umap'] = adata.obsm[umap_key_dif].copy()
                    
                    # Plot CellDiffusion - Batch (top left)
                    sc.pl.umap(
                        adata,
                        color=params.batch_key,
                        ax=axes[0, 0],
                        show=False,
                        frameon=False,
                        title=f"CellDiffusion - Batch"
                    )
                    
                    # Plot CellDiffusion - Labels (top right)
                    sc.pl.umap(
                        adata,
                        color=params.label_key,
                        ax=axes[0, 1],
                        show=False,
                        frameon=False,
                        title=f"CellDiffusion - Labels"
                    )
                    
                    # Clean up
                    del adata.obsm['X_umap']
                else:
                    print(f"    Warning: {umap_key_dif} not found in adata")
                    axes[0, 0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
                    axes[0, 1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            except Exception as e:
                print(f"    Error plotting CellDiffusion: {e}")
                import traceback
                traceback.print_exc()
                axes[0, 0].text(0.5, 0.5, 'Error loading data', ha='center', va='center')
                axes[0, 1].text(0.5, 0.5, 'Error loading data', ha='center', va='center')
        else:
            print(f"    Warning: CellDiffusion data not found for layer {layer}")
            axes[0, 0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            axes[0, 1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        
        # Plot GCN if available
        if 'gcn' in layer_data[layer]:
            try:
                umap_key_gcn = layer_data[layer]['gcn']
                print(f"    Plotting GCN using {umap_key_gcn}...")
                
                if umap_key_gcn in adata.obsm:
                    # Temporarily set UMAP for plotting
                    adata.obsm['X_umap'] = adata.obsm[umap_key_gcn].copy()
                    
                    # Plot GCN - Batch (bottom left)
                    sc.pl.umap(
                        adata,
                        color=params.batch_key,
                        ax=axes[1, 0],
                        show=False,
                        frameon=False,
                        title=f"GCN - Batch"
                    )
                    
                    # Plot GCN - Labels (bottom right)
                    sc.pl.umap(
                        adata,
                        color=params.label_key,
                        ax=axes[1, 1],
                        show=False,
                        frameon=False,
                        title=f"GCN - Labels"
                    )
                    
                    # Clean up
                    del adata.obsm['X_umap']
                else:
                    print(f"    Warning: {umap_key_gcn} not found in adata")
                    axes[1, 0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
                    axes[1, 1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            except Exception as e:
                print(f"    Error plotting GCN: {e}")
                import traceback
                traceback.print_exc()
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
