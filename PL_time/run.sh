#!/bin/bash

# Run script for PL_time timing benchmark pipeline
# Compares runtime of CellDiffusion, Harmony, and scVI integration methods

set -e

# Load miniconda module if available
if command -v module &> /dev/null; then
    echo "Loading miniconda module..."
    module load miniconda || true
fi

# Initialize conda
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
else
    echo "Error: conda not found. Please install miniconda/anaconda first."
    exit 1
fi

# Try to find and activate an appropriate environment
if conda env list | grep -qE "^\s*cdiff_time_env\s"; then
    echo "Activating conda environment: cdiff_time_env"
    conda activate cdiff_time_env
elif conda env list | grep -qE "^\s*cdiff_bbknn_env\s"; then
    echo "Activating conda environment: cdiff_bbknn_env"
    conda activate cdiff_bbknn_env
elif conda env list | grep -qE "^\s*cdiff_uniport_env\s"; then
    echo "Activating conda environment: cdiff_uniport_env"
    conda activate cdiff_uniport_env
elif conda env list | grep -qE "^\s*dif_snake_scib_env\s"; then
    echo "Activating conda environment: dif_snake_scib_env"
    conda activate dif_snake_scib_env
else
    echo "No suitable conda environment found. Available environments:"
    conda env list | grep -E "^\s+[a-zA-Z]"
    echo ""
    echo "Please create a suitable environment with CellDiffusion, Harmony, and scVI installed."
    exit 1
fi

# Run snakemake pipeline
echo ""
echo "=========================================="
echo "Running PL_time benchmark pipeline"
echo "=========================================="
echo ""

CORES=${1:-4}
echo "Using $CORES cores"

snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cores "$CORES" \
    --use-conda \
    --reason

echo ""
echo "=========================================="
echo "Pipeline complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Benchmark table: ../data/outputs/PL_time/timing_benchmark_table.csv"
echo "  Benchmark plot:  ../data/outputs/PL_time/timing_benchmark_plot.pdf"
echo ""
