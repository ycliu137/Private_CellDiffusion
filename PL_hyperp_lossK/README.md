# PL_hyperp_k: k_mnn Hyperparameter Optimization Pipeline

This Snakemake pipeline tests different `k_mnn` values for CellDiffusion integration and evaluates the results using scib-metrics.

## Overview

The pipeline:
1. Preprocesses single-cell data
2. Encodes features using an autoencoder
3. Builds integration graphs with different `k_mnn` values (10, 20, 30, 40, 50)
4. Runs CellDiffusion integration for each `k_mnn` value
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
- k_mnn values to test: `[10, 20, 30, 40, 50]`
- UMAP parameters
- Evaluation parameters

## Pipeline Structure

```
preprocess → encode_features → build_graph (for each k_mnn) → integrate (for each k_mnn)
                                                                  ↓
                                            aggregate_X_dif → compute_umap → plot_umap
                                                                  ↓
                                            scib_evaluation → scib_evaluation_plot
```

## Output Files

- `preprocessed.h5ad` - Preprocessed data
- `encoded.h5ad` - Feature encoded data (contains `X_fae`)
- `graph_kmnn{k_mnn}.h5ad` - Integration graphs for each k_mnn
- `integrated_kmnn{k_mnn}.h5ad` - Integrated data for each k_mnn (contains `X_dif`)
- `aggregated_X_dif.h5ad` - All X_dif embeddings aggregated (contains `X_dif_kmnn{k_mnn}` for each k_mnn)
- `data_with_umap.h5ad` - Aggregated data with UMAP embeddings (contains `X_umap_kmnn{k_mnn}`)
- `umap_plot.pdf` - UMAP visualizations for each k_mnn (colored by batch and labels)
- `scib_results_table.csv` - SCIB evaluation results table
- `scib_results_table_plot.pdf` - SCIB results table plot
- `scib_comparison_barplot.pdf` - Bar plot comparing three aggregate scores across k_mnn values

## Usage

```bash
cd PL_hyperp_k

# Edit config.yaml to set data paths and parameters
# Then run:
./run.sh

# Or run with custom options:
./run.sh --dry-run  # Preview what will be executed
./run.sh --jobs 4   # Use 4 cores instead of 8
```

## Key Parameters

- **k_mnn**: Number of mutual nearest neighbors for inter-batch connections in integration graph building
- **n_edges_per_node**: Maximum number of outgoing edges per node (default: 50)
- **Aggregate Scores**: Three main metrics compared across k_mnn values:
  - Total
  - Batch correction
  - Bio conservation

## Notes

- GPU resources are limited to prevent CUDA conflicts during parallel execution
- Each k_mnn value is processed independently, allowing parallel execution
- Results are aggregated for comparison and visualization
- UMAP is computed for each integrated embedding separately

