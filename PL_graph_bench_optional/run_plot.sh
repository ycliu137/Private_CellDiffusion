#!/usr/bin/env bash
set -euo pipefail

# run_plot.sh - Force re-run all plotting rules in the PL_graph_bench_optional pipeline
# Usage:
#   ./run_plot.sh            # execute plotting rules with 4 jobs (forced)
#   DRY_RUN=1 ./run_plot.sh  # perform snakemake dry-run
#   JOBS=8 ./run_plot.sh     # run with 8 cores

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

DRY_RUN=0
JOBS=4
FORCE=1

# Parse command-line options: -n (dry-run), -j N (jobs), -F (force=1), -f (force=0)
while getopts ":Fnj:f" opt; do
  case ${opt} in
    F ) FORCE=1 ;;
    n ) DRY_RUN=1 ;;
    j ) JOBS=${OPTARG} ;;
    f ) FORCE=0 ;;
    \? ) echo "Usage: $0 [-F] [-f] [-n] [-j jobs]"; exit 1 ;;
  esac
done
shift $((OPTIND -1))

# --- Environment setup (mirror run.sh) ---
if command -v module &> /dev/null; then
  module unload python3 2>/dev/null || true
  module load miniconda 2>/dev/null || true
fi

if command -v conda &> /dev/null; then
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

snakemake --unlock 2>/dev/null || true

# Plotting rules to run
RULES=(
  scib_evaluation
  scib_evaluation_plot
  scib_evaluation_plot_split
  plot_umap_combined
)

SNM_ARGS=( -s Snakefile --cores "$JOBS" )
if [[ "$DRY_RUN" == "1" ]]; then
  SNM_ARGS+=( -n )
fi

echo "Running plotting rules in PL_graph_bench_optional (DRY_RUN=${DRY_RUN}, JOBS=${JOBS}, FORCE=${FORCE})"

if [[ "$FORCE" == "1" ]]; then
  echo "Forcing execution of plotting rules (will run targets even if outputs exist)"
  snakemake "${SNM_ARGS[@]}" --forcerun "${RULES[@]}"
else
  snakemake "${SNM_ARGS[@]}" "${RULES[@]}"
fi

echo "Done."
