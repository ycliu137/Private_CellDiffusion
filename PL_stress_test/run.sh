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
    echo "Activating conda environment: dif_snake_scib_env"
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate dif_snake_scib_env 2>/dev/null || {
        echo "Warning: Could not activate dif_snake_scib_env. Continuing anyway..."
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
#   ./run.sh                    # Run with default settings (-F -j 8 -k --resources gpu=1)
#   ./run.sh --dry-run          # Dry run to see what will be executed
#   ./run.sh --jobs 4           # Run with up to 4 concurrent jobs
#   ./run.sh --jobs 1           # Run with 1 concurrent job (sequential)

# Default: use -F -j 8 -k with GPU resource limit
# -F: force re-execution of updated rules
# -j 8: use 8 cores for parallel execution
# -k: keep going even if some jobs fail
# --resources gpu=1: limit GPU access to prevent CUDA conflicts
if [ $# -eq 0 ]; then
    echo "Running snakemake with default settings: -F -j 8 -k --resources gpu=1"
    echo "Note: Using --resources gpu=1 to limit GPU access and prevent CUDA conflicts"
    snakemake -F -j 8 -k --resources gpu=1
else
    snakemake "$@"
fi

