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
    echo "Activating conda environment: dif_snake_scib_env"
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate dif_snake_scib_env 2>/dev/null || {
        echo "Error: Could not activate conda environment 'dif_snake_scib_env'"
        echo "Please create the environment first or update this script with the correct environment name"
        exit 1
    }
fi

# Unlock snakemake if needed
echo "Unlocking snakemake..."
snakemake --unlock 2>/dev/null || true

# Run snakemake with default or custom arguments
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

