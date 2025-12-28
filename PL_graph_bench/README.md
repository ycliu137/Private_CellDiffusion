# PL_graph_bench: Graph Building Method Benchmark Pipeline

This Snakemake pipeline benchmarks different graph building methods for CellDiffusion integration and evaluates the results using scib-metrics.

## Overview

The pipeline:
1. Preprocesses single-cell data
2. Encodes features using an autoencoder
3. Builds integration graphs using different methods (5 methods total)
4. Runs CellDiffusion integration for each graph building method
5. Aggregates all integration results
6. Evaluates integration quality using scib-metrics
7. Generates comparison bar plots for aggregate scores

## Graph Building Methods Tested

1. **build_integration_graph** - Default integration graph (uses Harmony embeddings)
2. **build_mnn_graph** - Simple MNN graph
3. **build_harmony_mnn_graph** - Harmony-based MNN graph
4. **build_scvi_mnn_graph** - SCVI-based MNN graph (trains SCVI model)
5. **build_omnn_scvi_graph** - Optimized SCVI-based MNN graph (trains SCVI model)

## Configuration

Edit `config.yaml` to set:
- Input data path
- Output directory
- Preprocessing parameters
- Feature encoder parameters
- Graph building methods to test (default: all 5 methods)
- Integration parameters (n_edges_per_node=50, k_mnn=50)
- Integration diffusion parameters
- Evaluation parameters
- SCIB evaluation parameters

## Pipeline Structure

```
preprocess → encode_features → build_graph (for each method) → integrate (for each method)
                                                                              ↓
                                            aggregate_X_dif → scib_evaluation → scib_evaluation_plot
```

## Output Files

- `preprocessed.h5ad` - Preprocessed data
- `encoded.h5ad` - Feature encoded data (contains `X_fae`)
- `graph_{method}.h5ad` - Integration graphs for each method
- `integrated_{method}.h5ad` - Integrated data for each method (contains `X_dif`)
- `aggregated_X_dif.h5ad` - All X_dif embeddings aggregated (contains `X_dif_{MethodName}` for each method)
- `scib_results_table.csv` - SCIB evaluation results table
- `scib_results_table_plot.pdf` - SCIB results table plot
- `scib_comparison_barplot.pdf` - Bar plot comparing three aggregate scores across graph building methods

## Usage

```bash
cd PL_graph_bench

# Edit config.yaml to set data paths and parameters
# Then run:
./run.sh

# Or run with custom options:
./run.sh --dry-run  # Preview what will be executed
./run.sh --jobs 4   # Use 4 cores instead of 8
```

## Key Parameters

- **Graph Building Methods**: Tested methods are defined in `config.yaml` under `graph_methods`
- **Common Parameters** (same for all methods):
  - `n_edges_per_node`: 50
  - `k_mnn`: 50
- **Aggregate Scores**: Three main metrics compared:
  - Total
  - Batch correction
  - Bio conservation

## Notes

- GPU resources are limited to prevent CUDA conflicts during parallel execution
- Each graph building method is processed independently, allowing parallel execution
- Some methods (SCVI-based) require training a model, which may take additional time
- Results are aggregated for comparison and visualization
- The bar plot shows performance comparison across all tested graph building methods

