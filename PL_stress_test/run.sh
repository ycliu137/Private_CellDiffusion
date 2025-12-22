#!/bin/bash
# Quick script to run the Snakemake pipeline

# Navigate to the pipeline directory
cd "$(dirname "$0")"

# Module management (for HPC environments)
if command -v module &> /dev/null; then
    echo "Unloading python3 module..."
    module unload python3 2>/dev/null || true
    
    echo "Loading miniconda module..."
    module load miniconda 2>/dev/null || true
fi

# Activate conda environment
if command -v conda &> /dev/null; then
    echo "Activating conda environment: dif_snake_env"
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate dif_snake_env 2>/dev/null || {
        echo "Warning: Could not activate dif_snake_env. Continuing anyway..."
    }
fi

# Unlock snakemake (in case of previous interruption)
echo "Unlocking snakemake..."
snakemake --unlock 2>/dev/null || true

# Check if snakemake is installed
if ! command -v snakemake &> /dev/null; then
    echo "Error: snakemake is not installed. Please install it first:"
    echo "  pip install snakemake"
    exit 1
fi

# Run snakemake with user-provided arguments or default
# Usage: ./run.sh [snakemake options]
# Examples:
#   ./run.sh                    # Run with default settings (--cores 8 --resources gpu=1)
#   ./run.sh --dry-run          # Dry run to see what will be executed
#   ./run.sh --cores 4          # Run with 4 cores
#   ./run.sh --cores 1          # Run with 1 core (sequential)

# Default: use 1 core with GPU resource limit to ensure GPU tasks run sequentially
# This prevents CUDA "device busy" errors when multiple tasks try to use GPU simultaneously
if [ $# -eq 0 ]; then
    echo "Running snakemake with default settings: --cores 1 --resources gpu=1"
    echo "Note: Using 1 core to ensure GPU tasks run sequentially and avoid CUDA conflicts"
    snakemake --cores 1 --resources gpu=1
else
    snakemake "$@"
fi

