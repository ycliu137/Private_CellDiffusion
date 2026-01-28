#!/bin/bash
# Build conda env for uniPort pipeline (integration + scib evaluation).
# Uses mamba. Optionally: module load miniconda (e.g. on HPC).
# Usage: ./uniport_env/build_env_mamba.sh

set -e

ENV_NAME="uniport_env"
PYTHON_VERSION="3.11"

echo "=========================================="
echo "Building uniPort environment with Mamba"
echo "=========================================="

if command -v module &> /dev/null; then
    echo "Loading miniconda module..."
    module load miniconda 2>/dev/null || true
fi

if ! command -v mamba &> /dev/null; then
    echo "Error: mamba not found. Install mamba or use: conda install mamba -n base -c conda-forge"
    exit 1
fi
echo "Mamba: $(which mamba)"

eval "$(conda shell.bash hook)"

if mamba env list | grep -qE "^\s*${ENV_NAME}\s"; then
    echo "Removing existing env: ${ENV_NAME}"
    mamba env remove -n ${ENV_NAME} -y
fi

echo "Creating ${ENV_NAME} with Python ${PYTHON_VERSION}"
mamba create -n ${ENV_NAME} python=${PYTHON_VERSION} -y

echo "Installing conda packages..."
mamba install -n ${ENV_NAME} -c conda-forge -y \
    numpy pandas scipy scikit-learn matplotlib networkx

echo "Installing PyTorch (GPU)..."
mamba install -n ${ENV_NAME} -c pytorch -c nvidia -c conda-forge -y \
    pytorch pytorch-cuda

echo "Installing pip packages (uniport, scanpy, scib, snakemake)..."
mamba run -n ${ENV_NAME} pip install \
    uniport \
    scanpy leidenalg python-igraph umap-learn \
    scvi-tools scib-metrics snakemake pyyaml ipykernel

echo ""
echo "=========================================="
echo "Environment build completed: ${ENV_NAME}"
echo "=========================================="
echo "  conda activate ${ENV_NAME}"
echo ""
