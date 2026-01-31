#!/bin/bash

# Build environment script for CellDiffusion + UniPort pipeline
# Uses mamba if available; falls back to conda.

set -e

ENV_NAME="cdiff_uniport_env"
PYTHON_VERSION="3.11"

echo "=========================================="
echo "Building CellDiffusion + UniPort Environment"
echo "=========================================="

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

# Choose package manager
if command -v mamba &> /dev/null; then
    PM="mamba"
    echo "Using mamba: $(which mamba)"
else
    PM="conda"
    echo "mamba not found, using conda: $(which conda)"
fi

# Remove existing environment if it exists
if $PM env list | grep -q "^${ENV_NAME}\s"; then
    echo "Removing existing environment: ${ENV_NAME}"
    $PM env remove -n ${ENV_NAME} -y
fi

# Create new environment
echo "Creating environment: ${ENV_NAME} with Python ${PYTHON_VERSION}"
$PM create -n ${ENV_NAME} python=${PYTHON_VERSION} -y

# Install core packages via conda-forge
echo "Installing core packages via conda-forge..."
$PM install -n ${ENV_NAME} -c conda-forge \
    numpy \
    pandas \
    "scipy<1.16" \
    scikit-learn \
    matplotlib \
    networkx \
    -y

# Install PyTorch (GPU version). If CPU-only, replace with pytorch-cpu.
echo "Installing PyTorch with GPU support..."
$PM install -n ${ENV_NAME} -c pytorch -c nvidia -c conda-forge \
    pytorch \
    pytorch-cuda \
    -y

# Install remaining packages via pip
echo "Installing pip dependencies..."
$PM run -n ${ENV_NAME} pip install -r requirements.txt
$PM run -n ${ENV_NAME} pip install -r uniport_env/requirements.txt

echo ""
echo "=========================================="
echo "Environment build completed!"
echo "=========================================="
echo ""
echo "To activate the environment, run:"
echo "  conda activate ${ENV_NAME}"
echo ""