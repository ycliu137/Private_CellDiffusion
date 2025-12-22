#!/bin/bash
# Quick script to run the Snakemake pipeline

# Navigate to the pipeline directory
cd "$(dirname "$0")"

# Check if snakemake is installed
if ! command -v snakemake &> /dev/null; then
    echo "Error: snakemake is not installed. Please install it first:"
    echo "  pip install snakemake"
    exit 1
fi

# Run snakemake with user-provided arguments or default
# Usage: ./run.sh [snakemake options]
# Examples:
#   ./run.sh                    # Run with default settings
#   ./run.sh --dry-run          # Dry run to see what will be executed
#   ./run.sh --cores 4          # Run with 4 cores
#   ./run.sh --cores 1 --rerun-incomplete  # Rerun incomplete jobs with 1 core

snakemake "$@"

