# PL_hyperp_OMNNK: n_edges_per_node Hyperparameter Optimization Pipeline

This Snakemake pipeline tests different `n_edges_per_node` values for CellDiffusion integration (with `k_mnn` fixed at 50) and evaluates the results using scib-metrics.

## Overview

The pipeline:
1. Preprocesses single-cell data
2. Encodes features using an autoencoder
3. Builds integration graphs with different `n_edges_per_node` values (10, 20, 30, 40, 50, 60, 70, 80, 90, 100), with `k_mnn` fixed at 50
4. Runs CellDiffusion integration for each `n_edges_per_node` value
5. Aggregates all integration results
6. Computes and visualizes UMAP embeddings
7. Evaluates integration quality using scib-metrics
8. Generates comparison bar plots for aggregate scores

## Configuration

Edit `config.yaml` to set:
- Input data path
- Output directory
- Preprocessing parameters
- Feature encoder parameters
- Integration parameters
- `n_edges_per_node` values to test: `[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]`
- `k_mnn` fixed value: `50`
- UMAP parameters
- Evaluation parameters

## Pipeline Structure

```
preprocess → encode_features → build_graph (for each n_edges_per_node) → integrate (for each n_edges_per_node)
                                                                              ↓
                                            aggregate_X_dif → compute_umap → plot_umap
                                                                              ↓
                                            scib_evaluation → scib_evaluation_plot
```

## Output Files

- `preprocessed.h5ad` - Preprocessed data
- `encoded.h5ad` - Feature encoded data (contains `X_fae`)
- `graph_nedges{n_edges_per_node}.h5ad` - Integration graphs for each n_edges_per_node
- `integrated_nedges{n_edges_per_node}.h5ad` - Integrated data for each n_edges_per_node (contains `X_dif`)
- `aggregated_X_dif.h5ad` - All X_dif embeddings aggregated (contains `X_dif_nedges{n_edges_per_node}` for each n_edges_per_node)
- `data_with_umap.h5ad` - Aggregated data with UMAP embeddings (contains `X_umap_nedges{n_edges_per_node}`)
- `umap_plot.pdf` - UMAP visualizations for each n_edges_per_node (colored by batch and labels)
- `scib_results_table.csv` - SCIB evaluation results table
- `scib_results_table_plot.pdf` - SCIB results table plot
- `scib_comparison_barplot.pdf` - Bar plot comparing three aggregate scores across n_edges_per_node values

## Usage

```bash
cd PL_hyperp_OMNNK

# Edit config.yaml to set data paths and parameters
# Then run:
./run.sh

# Or run with custom options:
./run.sh --dry-run  # Preview what will be executed
./run.sh --jobs 4   # Use 4 cores instead of 8
```

## Key Parameters

- **n_edges_per_node**: Maximum number of outgoing edges per node in the integration graph (varied: [10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
- **k_mnn**: Number of mutual nearest neighbors for inter-batch connections (fixed: 50)
- **Aggregate Scores**: Three main metrics compared across n_edges_per_node values:
  - Total
  - Batch correction
  - Bio conservation

## Notes

- GPU resources are limited to prevent CUDA conflicts during parallel execution
- Each n_edges_per_node value is processed independently, allowing parallel execution
- Results are aggregated for comparison and visualization
- UMAP is computed for each integrated embedding separately

