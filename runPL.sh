#!/bin/bash

# Run all PL_ pipelines sequentially
# Usage: ./runPL.sh [snakemake options to pass to each pipeline]
# Example: ./runPL.sh -j 8 --resources gpu=1
# 
# Note: This script continues running even if a pipeline fails,
#       and provides a summary at the end

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Unload conflicting Python module (if on HPC)
if command -v module &> /dev/null; then
    module unload python3 2>/dev/null || true
    module load miniconda 2>/dev/null || true
fi

# Activate conda environment (adjust name if needed)
if command -v conda &> /dev/null; then
    # Try to activate the environment
    if conda env list | grep -q "dif_snake_scib_env"; then
        echo "Activating conda environment: dif_snake_scib_env"
        conda activate dif_snake_scib_env
    elif conda env list | grep -q "dif_snake_env"; then
        echo "Activating conda environment: dif_snake_env"
        conda activate dif_snake_env
    else
        echo "Warning: Could not find conda environment. Please activate manually."
    fi
fi

# List of pipelines to run (in order)
PIPELINES=(
    "PL_hyperp"
    "PL_hyperp_k"
    "PL_hyperp_OMNNK"
    "PL_graph_bench"
    "PL_collapse_diag"
    "PL_stress_test"
    "PL_GCN"
)

# Function to run a single pipeline
run_pipeline() {
    local pipeline_dir="$1"
    shift  # Remove first argument, keep remaining as snakemake args
    
    if [ ! -d "$pipeline_dir" ]; then
        echo "Warning: Pipeline directory $pipeline_dir does not exist. Skipping."
        return 1
    fi
    
    if [ ! -f "$pipeline_dir/run.sh" ]; then
        echo "Warning: run.sh not found in $pipeline_dir. Skipping."
        return 1
    fi
    
    echo ""
    echo "=========================================="
    echo "Running pipeline: $pipeline_dir"
    echo "=========================================="
    echo ""
    
    # Change to pipeline directory and run
    cd "$SCRIPT_DIR/$pipeline_dir"
    
    # Run the pipeline's run.sh script with provided arguments
    if [ $# -eq 0 ]; then
        # No arguments, use default from run.sh
        ./run.sh
    else
        # Pass arguments to run.sh
        ./run.sh "$@"
    fi
    
    local exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        echo ""
        echo "✓ Pipeline $pipeline_dir completed successfully"
    else
        echo ""
        echo "✗ Pipeline $pipeline_dir failed with exit code $exit_code"
        echo "Continuing with remaining pipelines..."
    fi
    
    # Return to root directory
    cd "$SCRIPT_DIR"
    
    return $exit_code
}

# Main execution
echo "=========================================="
echo "Running All PL_ Pipelines"
echo "=========================================="
echo ""
echo "Pipelines to run:"
for pipeline in "${PIPELINES[@]}"; do
    echo "  - $pipeline"
done
echo ""

# Track results
SUCCESSFUL=0
FAILED=0
FAILED_PIPELINES=()

# Run each pipeline
for pipeline in "${PIPELINES[@]}"; do
    if run_pipeline "$pipeline" "$@"; then
        ((SUCCESSFUL++))
    else
        ((FAILED++))
        FAILED_PIPELINES+=("$pipeline")
    fi
    
    # Small delay between pipelines
    sleep 2
done

# Summary
echo ""
echo "=========================================="
echo "Summary"
echo "=========================================="
echo "Total pipelines: ${#PIPELINES[@]}"
echo "Successful: $SUCCESSFUL"
echo "Failed: $FAILED"

if [ $FAILED -gt 0 ]; then
    echo ""
    echo "Failed pipelines:"
    for pipeline in "${FAILED_PIPELINES[@]}"; do
        echo "  - $pipeline"
    done
    echo ""
    exit 1
else
    echo ""
    echo "All pipelines completed successfully!"
    exit 0
fi

