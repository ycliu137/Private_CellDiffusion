# CellDiffusion Stress Test Pipeline

This directory contains a Snakemake workflow for running the CellDiffusion stress test pipeline based on the `Data_Integration_stress_tests.ipynb` notebook.

## Directory Structure

```
PL_stress_test/
├── Snakefile              # Snakemake workflow definition
├── config.yaml            # Configuration file with all parameters
├── run.sh                 # Convenience script to run the pipeline
├── scripts/               # Python scripts for each pipeline step
│   ├── preprocess.py      # Data preprocessing
│   ├── encode_features.py # Feature encoding
│   └── integrate.py       # Integration pipeline
└── README.md              # This file
```

## Setup

1. Ensure you have the CellDiffusion environment activated:
   ```bash
   conda activate private_celldiffusion_env
   # or
   source activate private_celldiffusion_env
   ```

2. Install Snakemake if not already installed:
   ```bash
   pip install snakemake
   ```

3. Configure the paths and parameters in `config.yaml`:
   - Update `data.input_path` to point to your input h5ad file
   - Update `data.output_dir` to specify where results should be saved
   - Adjust any other parameters as needed

## Usage

### Quick Start

Use the provided run script:
```bash
cd PL_stress_test
./run.sh
```

Or run Snakemake directly:

### Basic usage

Run the complete pipeline:
```bash
cd PL_stress_test
snakemake
```

### Run with specific number of cores

```bash
snakemake --cores 1
# or using the run script:
./run.sh --cores 1
```

### Dry run (preview what will be executed)

```bash
snakemake --dry-run
# or:
./run.sh --dry-run
```

### Force re-run specific steps

```bash
snakemake --forceall integrate
```

### Clean intermediate files

```bash
snakemake --delete-all-output
```

### Run from project root

You can also run from the project root:
```bash
cd /path/to/Private_CellDiffusion
snakemake -s PL_stress_test/Snakefile
```

## Pipeline Steps

The pipeline consists of three main steps:

1. **Preprocessing** (`preprocess` rule)
   - Filter genes
   - Normalize and log-transform
   - Select highly variable genes

2. **Feature Encoding** (`encode_features` rule)
   - Train autoencoder to encode features
   - Creates `X_fae` representation

3. **Integration** (`integrate` rule)
   - Build integration loss adjacency
   - Build integration graph
   - Evaluate neighbor purity (before)
   - Add random edges
   - Run integration diffusion
   - Evaluate MNN neighbor purity (after)
   - Compute UMAP for visualization
   - Save results and metrics

## Output Files

- `data/outputs/stress_tests_pbmc10k/preprocessed.h5ad` - Preprocessed data
- `data/outputs/stress_tests_pbmc10k/encoded.h5ad` - Feature encoded data
- `data/outputs/stress_tests_pbmc10k/scDiffusion_integration_k{k}.h5ad` - Final integrated data for each k value (e.g., k=0, 10, 20, ..., 100)
- `data/outputs/stress_tests_pbmc10k/metrics_log.csv` - Metrics log with purity scores for all k values

## Configuration

All parameters can be adjusted in `config.yaml`. Key parameters include:

- **Preprocessing**: `min_cells`, `n_top_genes`, normalization parameters
- **Feature encoder**: Architecture dimensions, training epochs, learning rate
- **Integration graph**: Number of edges per node, MNN k parameter
- **Integration loss adjacency**: List of k values to test (e.g., `[0, 10, 20, ..., 100]`)
- **Random edges**: Number of random edges to add (`k_add`), random seed
- **Integration diffusion**: Model architecture, training parameters
- **Device**: `"cuda"` or `"cpu"`

## Notes

- The pipeline uses GPU if available and `device: "cuda"` is set in config
- Intermediate files are saved to allow resuming from any step
- Metrics are appended to the CSV file, so multiple runs will accumulate results
- The pipeline follows the exact workflow from `Data_Integration_stress_tests.ipynb`

