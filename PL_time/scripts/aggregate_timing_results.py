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

# Normalize dataset name (some scripts may store dataset as list)
def _normalize_dataset_name(value):
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return ""
        return value[0]
    return value

def _normalize_method_name(value):
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return ""
        return value[0]
    return value

def _normalize_numeric(value):
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return 0
        return value[0]
    return value

def _dataset_from_path(file_path):
    try:
        return Path(file_path).parent.name
    except Exception:
        return ""

# Create timing dataframe
timing_dfs = {}
for f, data in zip(timing_files, timing_data):
    dataset = _normalize_dataset_name(data['dataset'])
    if not dataset or dataset == "integration":
        dataset = _dataset_from_path(f)
    method = _normalize_method_name(data['method'])
    total_time = _normalize_numeric(data['total_time'])
    
    if dataset not in timing_dfs:
        timing_dfs[dataset] = {}
    
    timing_dfs[dataset][method] = {
        'total_time': total_time,
        'steps': data['steps']
    }

# Create stats dataframe
stats_dfs = {}
for f, data in zip(stats_files, stats_data):
    dataset = _normalize_dataset_name(data['dataset'])
    if not dataset or dataset == "integration":
        dataset = _dataset_from_path(f)
    method = _normalize_method_name(data['method'])
    
    if dataset not in stats_dfs:
        stats_dfs[dataset] = {}
    
    stats_dfs[dataset][method] = {
        'n_cells': _normalize_numeric(data['n_cells']),
        'n_batches': _normalize_numeric(data['n_batches'])
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
    
    # Add total time for each method (convert to minutes)
    for method in ['CellDiffusion', 'Harmony', 'scVI', 'Seurat']:
        if method in timing_info:
            total_time = _normalize_numeric(timing_info[method]['total_time'])
            # total_time is already in minutes in timing JSON
            row[f'{method}_Time(min)'] = round(total_time, 2)
        else:
            row[f'{method}_Time(min)'] = np.nan
    
    rows.append(row)

# Create dataframe
df = pd.DataFrame(rows)

# Ensure numeric columns are numeric (handle list/object types)
df['N_Cells'] = pd.to_numeric(df['N_Cells'], errors='coerce')
df['N_Batches'] = pd.to_numeric(df['N_Batches'], errors='coerce')

# Calculate speed ratios (relative to CellDiffusion)
if 'CellDiffusion_Time(min)' in df.columns:
    df['Harmony/CellDiff'] = round(df['Harmony_Time(min)'] / df['CellDiffusion_Time(min)'], 2)
    df['scVI/CellDiff'] = round(df['scVI_Time(min)'] / df['CellDiffusion_Time(min)'], 2)
    df['Seurat/CellDiff'] = round(df['Seurat_Time(min)'] / df['CellDiffusion_Time(min)'], 2)

# Add summary row
summary_row = {
    'Dataset': 'Average',
    'N_Cells': f"~{int(df['N_Cells'].mean(skipna=True))}",
    'N_Batches': f"~{df['N_Batches'].mean(skipna=True):.1f}",
}

for col in df.columns:
    if '_Time(min)' in col:
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

# Create visualization (table)
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    data_rows = df[df['Dataset'] != 'Average'].copy()
    if len(data_rows) > 0:
        table_cols = [
            'Dataset',
            'N_Cells',
            'N_Batches',
            'CellDiffusion_Time(min)',
            'Harmony_Time(min)',
            'scVI_Time(min)',
            'Seurat_Time(min)'
        ]

        # Format numeric columns for display
        display_df = data_rows[table_cols].copy()
        for col in ['N_Cells', 'N_Batches']:
            display_df[col] = display_df[col].astype(str)
        for col in ['CellDiffusion_Time(min)', 'Harmony_Time(min)', 'scVI_Time(min)', 'Seurat_Time(min)']:
            display_df[col] = display_df[col].map(lambda x: f"{x:.2f}" if pd.notnull(x) else "NA")

        n_rows = len(display_df)
        fig_height = max(2, 0.5 + 0.35 * n_rows)
        fig, ax = plt.subplots(figsize=(12, fig_height))
        ax.axis('off')

        table = ax.table(
            cellText=display_df.values,
            colLabels=display_df.columns,
            cellLoc='center',
            loc='center'
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.2)

        plt.tight_layout()
        Path(output_plot).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_plot, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Benchmark table plot saved to: {output_plot}")
except Exception as e:
    print(f"Warning: Could not create table plot: {e}")

print("\n=== Timing aggregation complete! ===")
