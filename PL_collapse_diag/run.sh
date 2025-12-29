#!/bin/bash

# Script to run the collapse diagnostic Snakemake pipeline
# This pipeline tests different network depths for CellDiffusion and GCN integration

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Unload python3 module if it exists (for HPC environments)
if command -v module &> /dev/null; then
    echo "Unloading python3 module..."
    module unload python3 2>/dev/null || true
    
    # Load miniconda module
    echo "Loading miniconda module..."
    module load miniconda 2>/dev/null || true
fi

# Activate conda environment
if [ -n "$CONDA_DEFAULT_ENV" ]; then
    echo "Conda environment already active: $CONDA_DEFAULT_ENV"
else
    echo "Activating conda environment: dif_snake_env"
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate dif_snake_env 2>/dev/null || {
        echo "Error: Could not activate conda environment 'dif_snake_env'"
        echo "Please create the environment first or update this script with the correct environment name"
        exit 1
    }
fi

# Unlock snakemake if needed
echo "Unlocking snakemake..."
snakemake --unlock 2>/dev/null || true

# Run snakemake with default or custom arguments
if [ $# -eq 0 ]; then
    echo "Running snakemake with default settings: -j 8 --resources gpu=1"
    echo "Note: Using --resources gpu=1 to limit GPU access and prevent CUDA conflicts"
    snakemake -j 8 --resources gpu=1
else
    echo "Running snakemake with custom arguments: $@"
    snakemake "$@"
fi

