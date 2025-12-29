#!/bin/bash

# Snakemake pipeline runner for k_mnn hyperparameter optimization
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

# Default: use -j 8 with GPU resource limit
# -j 8: use 8 cores for parallel execution
# --resources gpu=1: limit GPU access to prevent CUDA conflicts
# Note: -F (forceall) and -k (keep-going) are optional and can be added as arguments
# Usage examples:
#   ./run.sh                    # Run with default settings
#   ./run.sh -F -k              # Add -F (forceall) and -k (keep-going) to default settings
#   ./run.sh -j 4 -k            # Use 4 jobs with keep-going option
if [ $# -eq 0 ]; then
    echo "Running snakemake with default settings: -j 8 --resources gpu=1"
    echo "Note: Using --resources gpu=1 to limit GPU access and prevent CUDA conflicts"
    echo "Note: -F (forceall) and -k (keep-going) are optional. Use: ./run.sh -F -k"
    snakemake -j 8 --resources gpu=1
else
    echo "Running snakemake with custom arguments: $@"
    snakemake "$@"
fi

