# uniPort pipeline environment

Conda/pip environment for the PL_uniPort pipeline (uniPort integration + scib evaluation).

## Option 1: Conda from `environment.yml`

```bash
cd /path/to/Private_CellDiffusion
conda env create -f uniport_env/environment.yml
conda activate uniport_env
```

On CPU-only machines, edit `environment.yml` and remove or replace `pytorch-cuda` with `cpuonly`.

## Option 2: Mamba build script

```bash
./uniport_env/build_env_mamba.sh
conda activate uniport_env
```

On HPC with modules: `module load miniconda` before running the script if needed.

## Option 3: Pip into existing env

Install PyTorch (and CUDA if needed) first, then:

```bash
pip install -r uniport_env/requirements.txt
```

## Packages

- **uniport** – single-cell integration
- **scanpy**, **scvi-tools** – preprocessing and analysis
- **scib-metrics** – integration benchmarking
- **snakemake** – pipeline execution
