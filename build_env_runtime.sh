#!/bin/bash

# Build environment script for CellDiffusion + Harmony + scVI timing benchmark pipeline
# Uses mamba if available; falls back to conda.

set -e

ENV_NAME="cdiff_time_env"
PYTHON_VERSION="3.11"

echo "=========================================="
echo "Building PL_time Timing Benchmark Environment"
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
    anndata \
    umap-learn \
    -y

# Install PyTorch (GPU version). If CPU-only, replace with pytorch-cpu.
echo "Installing PyTorch with GPU support..."
$PM install -n ${ENV_NAME} -c pytorch -c nvidia -c conda-forge \
    pytorch \
    pytorch-cuda \
    -y

# Install scanpy, snakemake, and dependencies via conda-forge
echo "Installing scanpy, snakemake, and related packages..."
$PM install -n ${ENV_NAME} -c conda-forge \
    scanpy \
    snakemake \
    python-louvain \
    leidenalg \
    -y

# Install Harmony, scVI, and remaining packages via pip
echo "Installing Harmony, scVI, and additional dependencies..."
$PM run -n ${ENV_NAME} pip install harmonypy scvi-tools
$PM run -n ${ENV_NAME} pip install pytorch-lightning rich jax jaxlib numpyro
$PM run -n ${ENV_NAME} pip install -r requirements.txt

echo ""
echo "=========================================="
echo "Environment build completed!"
echo "=========================================="
echo ""
echo "To activate the environment, run:"
echo "  conda activate ${ENV_NAME}"
echo ""
echo "To run the timing benchmark pipeline:"
echo "  cd PL_time"
echo "  bash run.sh 4"
echo ""
