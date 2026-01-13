"""
Aggregate collapse diagnostic metrics from all evaluations
"""
import sys
from pathlib import Path
import pandas as pd

# Get input files from snakemake
input_csv_files = snakemake.input.csv_files
output_csv = snakemake.output.csv

print(f"\n=== Aggregating collapse metrics from {len(input_csv_files)} files ===")

# Read all CSV files and concatenate
all_results = []
for csv_file in input_csv_files:
    print(f"Loading: {csv_file}")
    df = pd.read_csv(csv_file)
    all_results.append(df)

# Combine all results
aggregated_df = pd.concat(all_results, ignore_index=True)

# Sort by method and network_layers
aggregated_df = aggregated_df.sort_values(['method', 'network_layers'])

print(f"\n=== Aggregated results ===")
print(aggregated_df.to_string(index=False))

# Save aggregated results
print(f"\n=== Saving aggregated results ===")
print(f"Saving to: {output_csv}")
Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
aggregated_df.to_csv(output_csv, index=False)

print("Aggregation complete!")

