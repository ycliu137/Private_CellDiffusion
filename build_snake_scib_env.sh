#!/bin/bash

# Build environment script for CellDiffusion using mamba
# This script uses module load miniconda and mamba to create the environment

set -e  # Exit on error

ENV_NAME="dif_snake_scib_env"
PYTHON_VERSION="3.12"

echo "=========================================="
echo "Building CellDiffusion Environment with Mamba"
echo "=========================================="

# Load miniconda module
echo ""
echo "Loading miniconda module..."
module load miniconda

# Check if mamba is available after loading module
if ! command -v mamba &> /dev/null; then
    echo "Error: mamba not found after loading miniconda module."
    exit 1
fi

echo "Mamba found: $(which mamba)"

# Initialize conda for bash shell (mamba uses conda's infrastructure)
eval "$(conda shell.bash hook)"

# Remove existing environment if it exists
echo ""
echo "Checking for existing environment..."
if mamba env list | grep -q "^${ENV_NAME}\s"; then
    echo "Removing existing environment: ${ENV_NAME}"
    mamba env remove -n ${ENV_NAME} -y
fi

# Create new environment with Python using mamba create
echo ""
echo "Creating new environment: ${ENV_NAME} with Python ${PYTHON_VERSION}"
mamba create -n ${ENV_NAME} python=${PYTHON_VERSION} -y

# Install core packages via conda-forge for better compatibility
echo ""
echo "Installing core packages via conda-forge..."
mamba install -n ${ENV_NAME} -c conda-forge \
    numpy=2.3.1 \
    pandas=2.3.1 \
    scipy=1.15.2 \
    scikit-learn=1.7.1 \
    matplotlib=3.10.3 \
    networkx=3.4.2 \
    -y

# Install PyTorch (GPU version with CUDA)
echo ""
echo "Installing PyTorch with GPU support..."
mamba install -n ${ENV_NAME} -c pytorch -c nvidia -c conda-forge \
    pytorch \
    pytorch-cuda \
    -y

# Install remaining packages via pip
echo ""
echo "Installing remaining packages via pip..."
mamba run -n ${ENV_NAME} pip install \
    leidenalg==0.10.2 \
    python-igraph==0.11.8 \
    python-louvain==0.16 \
    scanpy==1.11.3 \
    umap-learn==0.5.7 \
    scvi-tools \
    scib \
    ipykernel \
    snakemake
    
python -m ipykernel install --user --name ${ENV_NAME} --display-name "${ENV_NAME}"

echo ""
echo "=========================================="
echo "Environment build completed!"
echo "=========================================="
echo ""
echo "To activate the environment, run:"
echo "  module load miniconda"
echo "  conda activate ${ENV_NAME}"
echo ""

