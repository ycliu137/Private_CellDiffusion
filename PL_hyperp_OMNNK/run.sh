#!/bin/bash

# Snakemake pipeline runner for n_edges_per_node hyperparameter optimization
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

# Default: use -F -j 8 -k with GPU resource limit
# -F: force re-execution of updated rules
# -j 8: use 8 cores
# -k: keep going even if some jobs fail
if [ $# -eq 0 ]; then
    echo "Running snakemake with default settings: -F -j 8 -k --resources gpu=1"
    echo "Note: Using --resources gpu=1 to limit GPU access and prevent CUDA conflicts"
    snakemake -F -j 8 -k --resources gpu=1
else
    snakemake "$@"
fi

