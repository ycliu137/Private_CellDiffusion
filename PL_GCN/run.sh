#!/bin/bash

# Snakemake pipeline runner for CellDiffusion vs GCN comparison
# Usage: ./run.sh [snakemake options]

# Change to script directory
cd "$(dirname "$0")"

# Unload conflicting Python module (if on HPC)
if command -v module &> /dev/null; then
    module unload python3 2>/dev/null || true
    module load miniconda 2>/dev/null || true
fi

# Activate conda environment (adjust name if needed)
if command -v conda &> /dev/null; then
    # Try to activate the environment
    if conda env list | grep -q "dif_snake_scib_env"; then
        echo "Activating conda environment: dif_snake_scib_env"
        conda activate dif_snake_scib_env
    elif conda env list | grep -q "dif_snake_env"; then
        echo "Activating conda environment: dif_snake_env"
        conda activate dif_snake_env
    else
        echo "Warning: Could not find conda environment. Please activate manually."
    fi
fi

# Unlock snakemake in case of previous interruption
snakemake --unlock 2>/dev/null || true

# Default: use --jobs 1 with GPU resource limit to ensure GPU tasks run strictly sequentially
if [ $# -eq 0 ]; then
    echo "Running snakemake with default settings: --jobs 1 --resources gpu=1"
    echo "Note: Using --jobs 1 to ensure GPU tasks run strictly sequentially and avoid CUDA conflicts"
    snakemake --jobs 1 --resources gpu=1
else
    snakemake "$@"
fi

