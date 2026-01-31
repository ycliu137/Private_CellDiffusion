"""
Aggregate timing results from all methods and generate benchmark table.
Creates a comprehensive CSV with dataset info, timing for each method, and ratios.
"""
import sys
import json
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np

# Collect all timing and stats files
input_files = snakemake.input
output_table = snakemake.output.table
output_plot = snakemake.output.plot

print(f"\n=== Aggregating Timing Results ===")

# Group files by type
timing_files = [f for f in input_files if 'timing' in str(f)]
stats_files = [f for f in input_files if 'stats' in str(f)]

print(f"Found {len(timing_files)} timing files")
print(f"Found {len(stats_files)} stats files")

# Load all timing data
timing_data = []
for f in timing_files:
    with open(f, 'r') as file:
        data = json.load(file)
        timing_data.append(data)

# Load all stats data
stats_data = []
for f in stats_files:
    with open(f, 'r') as file:
        data = json.load(file)
        stats_data.append(data)

# Create timing dataframe
timing_dfs = {}
for data in timing_data:
    dataset = data['dataset']
    method = data['method']
    total_time = data['total_time']
    
    if dataset not in timing_dfs:
        timing_dfs[dataset] = {}
    
    timing_dfs[dataset][method] = {
        'total_time': total_time,
        'steps': data['steps']
    }

# Create stats dataframe
stats_dfs = {}
for data in stats_data:
    dataset = data['dataset']
    method = data['method']
    
    if dataset not in stats_dfs:
        stats_dfs[dataset] = {}
    
    stats_dfs[dataset][method] = {
        'n_cells': data['n_cells'],
        'n_batches': data['n_batches']
    }

# Build comprehensive benchmark table
rows = []
for dataset in sorted(timing_dfs.keys()):
    # Get stats from first method (should be same for all)
    stats_info = stats_dfs[dataset]
    first_method = list(stats_info.keys())[0]
    n_cells = stats_info[first_method]['n_cells']
    n_batches = stats_info[first_method]['n_batches']
    
    # Get timings for all methods
    timing_info = timing_dfs[dataset]
    
    row = {
        'Dataset': dataset,
        'N_Cells': n_cells,
        'N_Batches': n_batches,
    }
    
    # Add total time for each method
    for method in ['CellDiffusion', 'Harmony', 'scVI', 'Seurat']:
        if method in timing_info:
            row[f'{method}_Time(s)'] = round(timing_info[method]['total_time'], 2)
        else:
            row[f'{method}_Time(s)'] = np.nan
    
    rows.append(row)

# Create dataframe
df = pd.DataFrame(rows)

# Calculate speed ratios (relative to CellDiffusion)
if 'CellDiffusion_Time(s)' in df.columns:
    df['Harmony/CellDiff'] = round(df['Harmony_Time(s)'] / df['CellDiffusion_Time(s)'], 2)
    df['scVI/CellDiff'] = round(df['scVI_Time(s)'] / df['CellDiffusion_Time(s)'], 2)
    df['Seurat/CellDiff'] = round(df['Seurat_Time(s)'] / df['CellDiffusion_Time(s)'], 2)

# Add summary row
summary_row = {
    'Dataset': 'Average',
    'N_Cells': f"~{int(df['N_Cells'].mean())}",
    'N_Batches': f"~{df['N_Batches'].mean():.1f}",
}

for col in df.columns:
    if '_Time(s)' in col:
        summary_row[col] = round(df[col].mean(), 2)
    elif col.endswith('Diff'):
        summary_row[col] = round(df[col].mean(), 2)

df = pd.concat([df, pd.DataFrame([summary_row])], ignore_index=True)

# Save to CSV
Path(output_table).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output_table, index=False)
print(f"\nBenchmark table saved to: {output_table}")
print("\n" + "="*100)
print(df.to_string(index=False))
print("="*100)

# Create visualization
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    # Plot 1: Total time comparison
    data_rows = df[df['Dataset'] != 'Average'].copy()
    if len(data_rows) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Time comparison
        methods = ['CellDiffusion_Time(s)', 'Harmony_Time(s)', 'scVI_Time(s)', 'Seurat_Time(s)']
        method_labels = ['CellDiffusion', 'Harmony', 'scVI', 'Seurat']
        
        x = np.arange(len(data_rows))
        width = 0.2
        
        for i, method in enumerate(methods):
            if method in data_rows.columns:
                axes[0].bar(x + i*width, data_rows[method], width, label=method_labels[i])
        
        axes[0].set_xlabel('Dataset')
        axes[0].set_ylabel('Time (seconds)')
        axes[0].set_title('Integration Method Runtime Comparison')
        axes[0].set_xticks(x + width * 1.5)
        axes[0].set_xticklabels(data_rows['Dataset'], rotation=45, ha='right')
        axes[0].legend()
        axes[0].grid(axis='y', alpha=0.3)
        
        # Speed ratio
        if 'Harmony/CellDiff' in data_rows.columns and 'scVI/CellDiff' in data_rows.columns and 'Seurat/CellDiff' in data_rows.columns:
            axes[1].plot(data_rows['Dataset'], data_rows['Harmony/CellDiff'], 'o-', label='Harmony/CellDiff', linewidth=2, markersize=8)
            axes[1].plot(data_rows['Dataset'], data_rows['scVI/CellDiff'], 's-', label='scVI/CellDiff', linewidth=2, markersize=8)
            axes[1].plot(data_rows['Dataset'], data_rows['Seurat/CellDiff'], '^-', label='Seurat/CellDiff', linewidth=2, markersize=8)
            axes[1].axhline(y=1, color='k', linestyle='--', alpha=0.3, label='CellDiffusion baseline')
            axes[1].set_xlabel('Dataset')
            axes[1].set_ylabel('Relative Speed (ratio to CellDiffusion)')
            axes[1].set_title('Relative Speed Comparison')
            axes[1].set_xticklabels(data_rows['Dataset'], rotation=45, ha='right')
            axes[1].legend()
            axes[1].grid(alpha=0.3)
        
        plt.tight_layout()
        Path(output_plot).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_plot, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Benchmark plot saved to: {output_plot}")
except Exception as e:
    print(f"Warning: Could not create plot: {e}")

print("\n=== Timing aggregation complete! ===")
