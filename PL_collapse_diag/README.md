# PL_collapse_diag: Collapse Diagnostic Pipeline

This pipeline tests the effect of network depth on integration quality for both CellDiffusion and GCN integration methods.

## Overview

The pipeline:
1. Preprocesses single-cell data
2. Encodes features using an autoencoder
3. Builds integration graphs
4. Runs CellDiffusion integration with varying `num_steps_diffusion` = [2, 4, 6, 8, 10, 12, 14, 16]
5. Runs GCN integration with varying `num_layers_gcn` = [2, 4, 6, 8, 10, 12, 14, 16]
6. Evaluates neighbor purity for each configuration using `evaluate_knn_neighbor_purity`
7. Aggregates results and plots a line chart comparing both methods

## Pipeline Structure

```
PL_collapse_diag/
├── config.yaml           # Configuration file
├── Snakefile             # Snakemake workflow definition
├── run.sh                # Execution script
├── README.md             # This file
└── scripts/
    ├── preprocess.py              # Data preprocessing
    ├── encode_features.py         # Feature encoding
    ├── build_graph.py             # Graph building
    ├── integrate_celldiffusion.py # CellDiffusion integration (with num_steps wildcard)
    ├── integrate_gcn.py           # GCN integration (with num_layers wildcard)
    ├── evaluate_purity.py         # Neighbor purity evaluation
    ├── aggregate_metrics.py       # Aggregate all metrics
    └── plot_collapse_diag.py      # Plot line chart
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
- Graph data: `data_with_graph.h5ad`
- CellDiffusion results: `celldiffusion_nsteps{N}.h5ad` (for N in [2,4,6,8,10,12,14,16])
- GCN results: `gcn_nlayers{N}.h5ad` (for N in [2,4,6,8,10,12,14,16])
- Individual purity metrics: `purity_celldiffusion_nsteps{N}.csv` and `purity_gcn_nlayers{N}.csv`
- Aggregated metrics: `aggregated_metrics.csv`
- **Final plot**: `collapse_diagnostic_plot.pdf` - Line chart comparing neighbor purity vs network layers for both methods

## Plot Output

The final plot (`collapse_diagnostic_plot.pdf`) shows:
- X-axis: Network layers (2, 4, 6, 8, 10, 12, 14, 16)
- Y-axis: Neighbor purity (0.0 to 1.0)
- Two lines:
  - CellDiffusion (blue line with circles)
  - GCN (orange line with squares)

This visualization helps identify if and when network collapse occurs as network depth increases.

## Notes

- The pipeline uses GPU resources (`--resources gpu=1`) to limit concurrent GPU tasks
- Each integration run is independent, allowing for parallel execution
- Neighbor purity is evaluated using KNN graph built on the integrated embeddings
- The evaluation uses `evaluate_knn_neighbor_purity` from `sc_evaluation.neighbor_purity`

