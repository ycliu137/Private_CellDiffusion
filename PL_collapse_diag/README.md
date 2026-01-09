# PL_collapse_diag: Collapse Diagnostic Pipeline

This pipeline tests the effect of network depth on integration quality for both CellDiffusion and GCN integration methods.

## Overview

The pipeline:
1. Preprocesses single-cell data
2. Encodes features using an autoencoder
3. For each network layer:
   - Builds integration graph, runs CellDiffusion/GCN integration, and computes UMAP in one step (optimized to reduce file I/O)
   - CellDiffusion with varying `num_steps_diffusion` = [2, 4, 6, 8, 10, 12, 14, 16]
   - GCN with varying `num_layers_gcn` = [2, 4, 6, 8, 10, 12, 14, 16]
4. Plots UMAP visualizations (colored by batch and labels) for each configuration
8. Evaluates neighbor purity for each configuration using `evaluate_knn_neighbor_purity`
9. Aggregates all embeddings and runs SCIB evaluation
10. Generates multiple plots:
    - Neighbor purity vs network layers (line chart)
    - SCIB aggregate scores vs network layers (line chart with 3 subplots: Total, Batch correction, Bio conservation)

## Pipeline Structure

```
PL_collapse_diag/
├── config.yaml           # Configuration file
├── Snakefile             # Snakemake workflow definition
├── run.sh                # Execution script
├── README.md             # This file
└── scripts/
    ├── preprocess.py                      # Data preprocessing
    ├── encode_features.py                 # Feature encoding
    ├── build_integrate_umap_celldiffusion.py  # Build graph, integrate with CellDiffusion, and compute UMAP (optimized)
    ├── build_integrate_umap_gcn.py        # Build graph, integrate with GCN, and compute UMAP (optimized)
    ├── plot_umap_combined.py              # Plot combined UMAP (multi-page PDF)
    ├── evaluate_purity.py                 # Neighbor purity evaluation
    ├── evaluate_purity_mesenchymal.py     # Neighbor purity evaluation (Mesenchymal lineage only)
    ├── aggregate_metrics.py               # Aggregate purity metrics
    ├── aggregate_metrics_mesenchymal.py   # Aggregate purity metrics (Mesenchymal lineage only)
    ├── aggregate_embeddings.py            # Aggregate all embeddings for SCIB
    ├── scib_evaluation.py                 # SCIB benchmark evaluation
    ├── scib_evaluation_mesenchymal.py     # SCIB benchmark evaluation (Mesenchymal lineage only)
    ├── scib_evaluation_plot.py            # Plot SCIB aggregate scores line chart
    ├── scib_evaluation_plot_mesenchymal.py # Plot SCIB aggregate scores line chart (Mesenchymal lineage only)
    ├── plot_collapse_diag.py              # Plot neighbor purity line chart
    └── plot_collapse_diag_mesenchymal.py  # Plot neighbor purity line chart (Mesenchymal lineage only)
```

## Configuration

Edit `config.yaml` to set:
- Data paths (input and output directories)
- Preprocessing parameters
- Feature encoder parameters
- Integration parameters
- Evaluation parameters (label_key, batch_key, k for KNN)

## Usage

### Basic Usage

```bash
cd PL_collapse_diag
./run.sh
```

### Custom Snakemake Options

```bash
# Run with custom options
./run.sh -j 4 --resources gpu=1 -k

# Run specific targets
./run.sh collapse_diagnostic_plot.pdf
```

## Output Files

The pipeline generates:
- Preprocessed data: `preprocessed.h5ad`
- Encoded features: `encoded.h5ad`
- CellDiffusion results: `celldiffusion_nsteps{N}.h5ad` (for N in [2,4,6,8,10,12,14,16]) - contains graph, X_dif, and UMAP
- GCN results: `gcn_nlayers{N}.h5ad` (for N in [2,4,6,8,10,12,14,16]) - contains graph, X_gcn, and UMAP
- **UMAP plots**: `umap_combined_multipage.pdf` - Multi-page PDF with combined CellDiffusion and GCN UMAPs for each network layer
- Individual purity metrics: `purity_celldiffusion_nsteps{N}.csv` and `purity_gcn_nlayers{N}.csv`
- Aggregated metrics: `aggregated_metrics.csv`
- Aggregated embeddings: `aggregated_embeddings.h5ad` (contains all X_dif and X_gcn embeddings)
- SCIB results: `scib_results_table.csv`, `scib_benchmarker.pkl`, `scib_results_table_plot.pdf`
- **Final plots**:
  - `collapse_diagnostic_plot.pdf` - Neighbor purity vs network layers (line chart)
  - `scib_evaluation_lineplot.pdf` - SCIB aggregate scores vs network layers (3 subplots: Total, Batch correction, Bio conservation)

## Plot Outputs

### 1. Neighbor Purity Plot (`collapse_diagnostic_plot.pdf`)
- X-axis: Network layers (2, 4, 6, 8, 10, 12, 14, 16)
- Y-axis: Neighbor purity (0.0 to 1.0)
- Two lines:
  - CellDiffusion (blue line with circles)
  - GCN (orange line with squares)

### 2. SCIB Aggregate Scores Plot (`scib_evaluation_lineplot.pdf`)
- Three subplots showing:
  - **Total**: Overall integration quality score
  - **Batch correction**: Batch correction quality score
  - **Bio conservation**: Biological information preservation score
- Each subplot:
  - X-axis: Network layers (2, 4, 6, 8, 10, 12, 14, 16)
  - Y-axis: Score (0.0 to 1.0)
  - Two lines: CellDiffusion (blue) and GCN (orange)

### 3. UMAP Plots
- **Combined multi-page PDF**: `umap_combined_multipage.pdf`
  - Each page corresponds to one network layer
  - Each page contains 4 subplots:
    - Top row: CellDiffusion (Batch, Labels)
    - Bottom row: GCN (Batch, Labels)
  - Allows easy comparison between methods for the same network depth
  - Reads UMAP data directly from `aggregated_embeddings.h5ad` to reduce file I/O

These visualizations help identify if and when network collapse occurs as network depth increases.

## Notes

- The pipeline uses GPU resources (`--resources gpu=1`) to limit concurrent GPU tasks
- Each integration run is independent, allowing for parallel execution
- Neighbor purity is evaluated using KNN graph built on the integrated embeddings
- The evaluation uses `evaluate_knn_neighbor_purity` from `sc_evaluation.neighbor_purity`

