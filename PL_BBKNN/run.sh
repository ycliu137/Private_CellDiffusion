#!/bin/bash
# Run PL_BBKNN pipeline (BBKNN integration + scib evaluation).
# Usage: ./run.sh [snakemake options]
# Requires cdiff_bbknn_env. Create via: ../build_env_bbknn.sh

cd "$(dirname "$0")"

if command -v module &> /dev/null; then
    module unload python3 2>/dev/null || true
    module load miniconda 2>/dev/null || true
fi

if command -v conda &> /dev/null; then
    if conda env list | grep -qE "^\s*cdiff_bbknn_env\s"; then
        echo "Activating conda environment: cdiff_bbknn_env"
        conda activate cdiff_bbknn_env
    elif conda env list | grep -qE "^\s*cdiff_uniport_env\s"; then
        echo "Activating conda environment: cdiff_uniport_env"
        conda activate cdiff_uniport_env
    elif conda env list | grep -q "dif_snake_scib_env"; then
        echo "Activating conda environment: dif_snake_scib_env"
        conda activate dif_snake_scib_env
    else
        echo "Warning: cdiff_bbknn_env not found. Create with: ../build_env_bbknn.sh"
    fi
fi

snakemake --unlock 2>/dev/null || true

if [ $# -eq 0 ]; then
    echo "Running snakemake: -j 4"
    snakemake -j 4
else
    echo "Running snakemake: $*"
    snakemake "$@"
fi
