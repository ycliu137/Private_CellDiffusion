# PL_uniPort

Pipeline for **uniPort** integration and **scib** evaluation of scRNA-seq datasets.

## Datasets

- Works with H5AD format scRNA-seq data files
- Each dataset must be a single H5AD file containing multiple batches/datasets
- Batches are identified via the `batch_key` column in `adata.obs`

## Setup

1. **Create environment** (from project root):

   ```bash
   # Option A: mamba
   ./build_env_mamba.sh

   # Option B: conda  
   conda env create -f uniport_env/environment.yml
   ```

   Then activate:

   ```bash
   conda activate uniport_env
   ```

2. **Configure** `config.yaml`:

   - `datasets`: List of dataset names (H5AD files in `input_path/{dataset}.h5ad`)
   - `data.input_path`: Directory containing H5AD files
   - `data.output_dir`: Output directory for results
   - `evaluation.batch_key`: Column name for batch labels (default: "batch")
   - `evaluation.label_key`: Column name for cell type labels (default: "labels")
   - `uniport.*`: uniPort settings (HVG counts, mode, lambda_s, etc.)

## Run

```bash
cd PL_uniPort
./run.sh
```

Or run Snakemake directly:

```bash
conda activate uniport_env
cd PL_uniPort
snakemake -j 4
```

## Outputs

Per-dataset outputs in `{output_dir}/{dataset}/`:

- `uniport_integrated.h5ad` – integrated AnnData with `obsm['X_uniport']` (latent representation)
- `aggregated_embeddings.h5ad` – same as above (for pipeline consistency)
- `scib_benchmarker.pkl` – SCIB benchmarker object
- `scib_results_table.csv` – metrics table (ARI, NMI, ASW, silhouette, kBET, LISI, etc.)
- `scib_results_table_plot.pdf` – visualization of SCIB metrics

## Example Configuration

```yaml
datasets: ["dataset1", "dataset2", "dataset3"]

data:
  input_path: "../data/inputs/integration/"
  output_dir: "../data/outputs/PL_uniPort"

uniport:
  n_hvg_common: 2000      # Common HVGs
  n_hvg_specific: 2000    # Per-batch specific HVGs
  mode: "h"               # UniPort mode
  lambda_s: 1.0           # Scaling parameter

evaluation:
  batch_key: "batch"      # Column with batch/dataset labels
  label_key: "labels"     # Column with cell type labels
```

## Environment files

Env assets live in **`uniport_env/`**:

- `environment.yml` – conda env spec.
- `requirements.txt` – pip-only deps.
- `build_env_mamba.sh` – build script.

