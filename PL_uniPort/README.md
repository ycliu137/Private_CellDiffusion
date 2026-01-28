# PL_uniPort

Pipeline for **uniPort** integration and **scib** evaluation, based on `Uniport_integration.py`.

## Setup

1. **Create environment** (from project root):

   ```bash
   # Option A: mamba
   ./uniport_env/build_env_mamba.sh

   # Option B: conda
   conda env create -f uniport_env/environment.yml
   ```

   Then activate:

   ```bash
   conda activate uniport_env
   ```

2. **Configure** `config.yaml`:

   - `data.datasets`: list of `{path, prefix}` for 10x MTX folders (each `path` contains `matrix.mtx`, `genes.tsv`, `barcodes.tsv`).
   - `data.output_dir`: where to write outputs.
   - `uniport.*`: uniPort settings (HVG, `mode`, `lambda_s`, etc.).
   - `evaluation.batch_key` / `label_key`: `adata.obs` keys for scib (default `source`, `domain_id`).

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

- `{output_dir}/uniport_integrated.h5ad` – integrated AnnData (`obsm['X_uniport']` = latent).
- `{output_dir}/figs/` – UMAP overview plots.
- `{output_dir}/scib_benchmarker.pkl` – scib benchmarker object.
- `{output_dir}/scib_results_table.csv` – metrics table.
- `{output_dir}/scib_results_table_plot.pdf` – scib heatmap.
- `{output_dir}/scib_comparison_barplot.pdf` – bar plot of aggregate scores.

## Environment files

Env assets live in **`uniport_env/`**:

- `environment.yml` – conda env spec.
- `requirements.txt` – pip-only deps.
- `build_env_mamba.sh` – build script.
- `README.md` – usage notes.
