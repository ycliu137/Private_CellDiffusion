#!/usr/bin/env bash
set -euo pipefail

# run_plot.sh - Run all plotting rules in the PL_collapse_diag Snakemake pipeline
# Usage:
#   ./run_plot.sh            # execute plotting rules with 4 jobs
#   DRY_RUN=1 ./run_plot.sh  # perform snakemake dry-run
#   JOBS=8 ./run_plot.sh     # run with 8 cores

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

DRY_RUN=0
JOBS=4
FORCE=0

# Parse command-line options: -F (force), -n (dry-run), -j N (jobs)
while getopts ":Fnj:" opt; do
  case ${opt} in
    F ) FORCE=1 ;;
    n ) DRY_RUN=1 ;;
    j ) JOBS=${OPTARG} ;;
    \? ) echo "Usage: $0 [-F] [-n] [-j jobs]"; exit 1 ;;
  esac
done
shift $((OPTIND -1))

# List of plotting rules to run
RULES=(
  plot_collapse_diag
  plot_collapse_metrics
  plot_umap_network_layers_comparison
  plot_collapse_diag_mesenchymal
  scib_evaluation_plot
  scib_evaluation_plot_mesenchymal
)

# Build snakemake target arguments
TARGETS=()
for r in "${RULES[@]}"; do
  TARGETS+=("$r")
done

SNM_ARGS=( -s Snakefile --cores "$JOBS" )
if [[ "$DRY_RUN" == "1" ]]; then
  SNM_ARGS+=( -n )
fi

# Report execution mode
echo "Running plotting rules in PL_collapse_diag (DRY_RUN=${DRY_RUN}, JOBS=${JOBS}, FORCE=${FORCE})"
# Run snakemake with the selected targets
if [[ "$FORCE" == "1" ]]; then
  echo "Forcing execution of plotting rules (will run targets even if outputs exist)"
  snakemake "${SNM_ARGS[@]}" --forcerun "${RULES[@]}"
else
  snakemake "${SNM_ARGS[@]}" "${TARGETS[@]}"
fi

echo "Done."
