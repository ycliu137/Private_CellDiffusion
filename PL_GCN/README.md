# PL_GCN: CellDiffusion vs GCN Comparison Pipeline

This Snakemake pipeline compares CellDiffusion and GCN integration methods on single-cell RNA-seq data.

## Overview

The pipeline:
1. Preprocesses single-cell data
2. Trains a feature encoder (autoencoder) to generate `X_fae` embeddings
3. Builds an integration graph and loss adjacency
4. Runs CellDiffusion integration → saves embeddings as `X_dif`
5. Runs GCN integration → saves embeddings as `X_gcn`
6. Computes UMAP for both embeddings → saves as `X_umap_dif` and `X_umap_gcn`
7. Plots UMAP visualizations
8. Benchmarks both methods using `scib-metrics`

## Directory Structure

```
PL_GCN/
├── config.yaml          # Configuration file
├── Snakefile           # Snakemake workflow definition
├── run.sh              # Pipeline runner script
├── README.md           # This file
└── scripts/
    ├── preprocess.py
    ├── encode_features.py
    ├── build_graph.py
    ├── integrate_celldiffusion.py
    ├── integrate_gcn.py
    ├── compute_umap.py
    ├── plot_umap.py
    └── scib_evaluation.py
```

## Configuration

Edit `config.yaml` to set:
- Input data path
- Output directory
- Preprocessing parameters
- Feature encoder parameters
- Integration parameters (for both CellDiffusion and GCN)
- UMAP parameters
- Evaluation parameters
- Device (cuda/cpu)

## Usage

### Basic usage

```bash
cd PL_GCN
./run.sh
```

### Custom options

```bash
# Dry run (see what will be executed)
./run.sh --dry-run

# Run with specific number of jobs
./run.sh --jobs 2

# Force rerun specific rules
./run.sh --forceall
```

## Output Files

- `preprocessed.h5ad` - Preprocessed data
- `encoded.h5ad` - Data with `X_fae` embeddings
- `data_with_graph.h5ad` - Data with integration graph
- `celldiffusion_integrated.h5ad` - Data with `X_dif` embeddings
- `gcn_integrated.h5ad` - Data with `X_gcn` embeddings
- `final_integrated.h5ad` - Combined data with both `X_dif` and `X_gcn`
- `data_with_umap.h5ad` - Data with UMAP embeddings
- `umap_plot.pdf` - UMAP visualizations
- `scib_results_table.csv` - SCIB benchmark results
- `scib_benchmarker.pkl` - Full benchmarker object
- `scib_results_table_plot.pdf` - SCIB results visualization

## Requirements

- Python 3.11+
- Snakemake
- Conda/Mamba environment with:
  - scanpy
  - celldiffusion
  - scib-metrics
  - pytorch (with CUDA support if using GPU)
  - matplotlib
  - pandas
  - numpy

## Notes

- The pipeline uses GPU resources with `resources: gpu=1` to prevent conflicts
- Default execution uses `--jobs 1` to ensure sequential GPU execution
- All embeddings are saved in `adata.obsm`:
  - `X_dif`: CellDiffusion embeddings
  - `X_gcn`: GCN embeddings
  - `X_umap_dif`: UMAP of CellDiffusion embeddings
  - `X_umap_gcn`: UMAP of GCN embeddings

## References

Based on:
- `Data_Integration.ipynb` - CellDiffusion integration workflow
- `GCN_Integration.ipynb` - GCN integration workflow

