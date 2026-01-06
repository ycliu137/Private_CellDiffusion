# PL_hyperp_lossK: integration_loss_adj.k Hyperparameter Optimization Pipeline

This Snakemake pipeline tests different `integration_loss_adj.k` (loss_k) values for CellDiffusion integration and evaluates the results using scib-metrics.

## Overview

The pipeline:
1. Preprocesses single-cell data
2. Encodes features using an autoencoder
3. Builds integration graphs with a fixed `k_mnn` value (default: 50) but different `loss_k` values for loss adjacency
4. Runs CellDiffusion integration for each `loss_k` value
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
- `k_mnn`: Fixed value for integration graph building (default: 50)
- `integration_loss_adj.k`: List of loss_k values to test (e.g., `[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]`)
- UMAP parameters
- Evaluation parameters

## Pipeline Structure

```
preprocess → encode_features → build_graph (for each loss_k) → integrate (for each loss_k)
                                                                  ↓
                                            aggregate_X_dif → compute_umap → plot_umap
                                                                  ↓
                                            scib_evaluation → scib_evaluation_plot
```

## Output Files

- `preprocessed.h5ad` - Preprocessed data
- `encoded.h5ad` - Feature encoded data (contains `X_fae`)
- `graph_lossk{loss_k}.h5ad` - Integration graphs for each loss_k
- `integrated_lossk{loss_k}.h5ad` - Integrated data for each loss_k (contains `X_dif`)
- `aggregated_X_dif.h5ad` - All X_dif embeddings aggregated (contains `X_dif_lossk{loss_k}` for each loss_k)
- `data_with_umap.h5ad` - Aggregated data with UMAP embeddings (contains `X_umap_lossk{loss_k}`)
- `umap_plot.pdf` - UMAP visualizations for each loss_k (colored by batch and labels)
- `scib_results_table.csv` - SCIB evaluation results table
- `scib_results_table_plot.pdf` - SCIB results table plot
- `scib_comparison_barplot.pdf` - Bar plot comparing three aggregate scores across loss_k values

## Usage

```bash
cd PL_hyperp_lossK

# Edit config.yaml to set data paths and parameters
# Then run:
./run.sh

# Or run with custom options:
./run.sh --dry-run  # Preview what will be executed
./run.sh --jobs 4   # Use 4 cores instead of 8
```

## Key Parameters

- **k_mnn**: Fixed number of mutual nearest neighbors for inter-batch connections in integration graph building (default: 50)
- **integration_loss_adj.k (loss_k)**: List of k values for loss adjacency in CellDiffusion integration (varied in this pipeline)
- **n_edges_per_node**: Maximum number of outgoing edges per node (default: 50)
- **Aggregate Scores**: Three main metrics compared across loss_k values:
  - Total
  - Batch correction
  - Bio conservation

## Notes

- GPU resources are limited to prevent CUDA conflicts during parallel execution
- Each loss_k value is processed independently, allowing parallel execution
- Results are aggregated for comparison and visualization
- UMAP is computed for each integrated embedding separately

