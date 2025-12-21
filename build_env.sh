#!/bin/bash

# Build environment script for CellDiffusion
# This script creates a new conda/mamba environment for the CellDiffusion project

set -e  # Exit on error

ENV_NAME="CellDiffusion_env"
PYTHON_VERSION="3.11"

echo "=========================================="
echo "Building CellDiffusion Environment"
echo "=========================================="

# Check if micromamba is available, otherwise use mamba or conda
if command -v micromamba &> /dev/null; then
    CONDA_CMD="micromamba"
    echo "Using micromamba"
elif command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "Using mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    echo "Using conda"
else
    echo "Error: Neither micromamba, mamba, nor conda found. Please install one of them."
    exit 1
fi

# Remove existing environment if it exists
echo ""
echo "Checking for existing environment..."
if $CONDA_CMD env list | grep -q "^${ENV_NAME}\s"; then
    echo "Removing existing environment: ${ENV_NAME}"
    $CONDA_CMD env remove -n ${ENV_NAME} -y
fi

# Create new environment with Python
echo ""
echo "Creating new environment: ${ENV_NAME} with Python ${PYTHON_VERSION}"
$CONDA_CMD create -n ${ENV_NAME} python=${PYTHON_VERSION} -y

# Activate environment (use source for bash/zsh compatibility)
echo ""
echo "Activating environment..."
if [ "$CONDA_CMD" = "micromamba" ]; then
    eval "$($CONDA_CMD shell hook --shell bash)"
    $CONDA_CMD activate ${ENV_NAME}
else
    source $($CONDA_CMD info --base)/etc/profile.d/conda.sh
    conda activate ${ENV_NAME}
fi

# Install core packages via conda-forge for better compatibility
echo ""
echo "Installing core packages via conda-forge..."
$CONDA_CMD install -n ${ENV_NAME} -c conda-forge \
    numpy=2.3.1 \
    pandas=2.3.1 \
    scipy=1.15.2 \
    scikit-learn=1.7.1 \
    matplotlib=3.10.3 \
    networkx=3.4.2 \
    -y

# Install PyTorch (CPU version)
echo ""
echo "Installing PyTorch..."
$CONDA_CMD install -n ${ENV_NAME} -c pytorch -c conda-forge pytorch cpuonly -y

# Install remaining packages via pip
echo ""
echo "Installing remaining packages via pip..."
$CONDA_CMD run -n ${ENV_NAME} pip install \
    leidenalg==0.10.2 \
    python-igraph==0.11.8 \
    python-louvain==0.16 \
    scanpy==1.11.3 \
    umap-learn==0.5.7 \
    scvi-tools \
    ipykernel

echo ""
echo "=========================================="
echo "Environment build completed!"
echo "=========================================="
echo ""
echo "To activate the environment, run:"
if [ "$CONDA_CMD" = "micromamba" ]; then
    echo "  micromamba activate ${ENV_NAME}"
else
    echo "  conda activate ${ENV_NAME}"
fi
echo ""

