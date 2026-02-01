"""
Plot SCIB evaluation results: separate plots for Total, Batch correction, Bio conservation.
Each plot compares CellDiffusion vs GCN across graph building methods.
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Global plot style for publication-quality figures
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "axes.linewidth": 1.0,
    "figure.dpi": 300,
})

# Get input and output from snakemake
input_table = snakemake.input.table
output_total = snakemake.output.total_pdf
output_batch = snakemake.output.batch_pdf
output_bio = snakemake.output.bio_pdf
output_combined = snakemake.output.combined_pdf

print("\n=== Loading SCIB results ===")
print(f"Input file: {input_table}")

# Read CSV file
df = pd.read_csv(input_table, index_col=0)
print(f"Results table shape: {df.shape}")
print(f"Columns: {list(df.columns)}")
print(f"\nIndex (Embedding methods): {list(df.index)}")

# Check for "Metric Type" row and remove it if present
METRIC_TYPE_ROW = "Metric Type"
if METRIC_TYPE_ROW in df.index:
    print(f"\nRemoving '{METRIC_TYPE_ROW}' row from data")
    df = df.drop(METRIC_TYPE_ROW)

# Find all graph method rows (X_dif_{method} and X_gcn_{method})
print("\n=== Finding integration methods in results ===")
method_rows = {}  # {(method, integration_type): row_name}
for method_name in df.index:
    method_str = str(method_name)
    if method_str.startswith('X_dif_'):
        method_key = method_str.replace('X_dif_', '')
        method_rows[(method_key, 'CellDiffusion')] = method_name
        print(f"Found CellDiffusion method: {method_key} -> {method_name}")
    elif method_str.startswith('X_gcn_'):
        method_key = method_str.replace('X_gcn_', '')
        method_rows[(method_key, 'GCN')] = method_name
        print(f"Found GCN method: {method_key} -> {method_name}")

if len(method_rows) == 0:
    raise ValueError("Could not find any X_dif or X_gcn methods in results table. Available methods: " + str(list(df.index)))

# Extract unique graph building methods and integration types
all_method_keys = sorted(set(method for method, _ in method_rows.keys()))
integration_types = ['CellDiffusion', 'GCN']
print(f"\nFound {len(method_rows)} integration methods across {len(all_method_keys)} graph building methods")
print(f"Graph building methods: {all_method_keys}")
print(f"Integration types: {integration_types}")

# Extract the three aggregate score columns
aggregate_score_cols = {
    'Total': None,
    'Batch correction': None,
    'Bio conservation': None
}

print("\n=== Extracting aggregate scores ===")
for col in df.columns:
    col_lower = col.lower()
    if col == 'Total' or 'total' in col_lower:
        aggregate_score_cols['Total'] = col
    elif 'batch' in col_lower and 'correction' in col_lower:
        aggregate_score_cols['Batch correction'] = col
    elif ('bio' in col_lower or 'conservation' in col_lower) and 'bio' in col_lower:
        aggregate_score_cols['Bio conservation'] = col

found_scores = {k: v for k, v in aggregate_score_cols.items() if v is not None}
print("Found aggregate score columns:")
for score_name, col_name in found_scores.items():
    print(f"  {score_name}: {col_name}")

missing_scores = [k for k, v in aggregate_score_cols.items() if v is None]
if missing_scores:
    raise ValueError("Missing aggregate score columns: " + str(missing_scores) + ". Available columns: " + str(list(df.columns)))

# Convert numeric columns to float (in case they were read as strings)
for col_name in found_scores.values():
    df[col_name] = pd.to_numeric(df[col_name], errors='coerce')

# Extract scores for each method and integration type
scores_data = {}  # {(method, integration_type): {score_name: value}}
for (method_key, integration_type), method_name in method_rows.items():
    key = (method_key, integration_type)
    scores_data[key] = {}
    for score_name, col_name in found_scores.items():
        value = df.loc[method_name, col_name]
        if pd.isna(value):
            print(f"  Warning: {score_name} is NaN for {method_key} ({integration_type})")
            value = 0.0
        scores_data[key][score_name] = float(value)
        print(f"  {method_key} ({integration_type}) - {score_name}: {value:.4f}")


def plot_single_score(score_name: str, output_pdf: str, colors: tuple[str, str]) -> None:
    print(f"\n=== Creating plot for {score_name} ===")

    # Set up plot
    fig, ax = plt.subplots(figsize=(8.5, 5))

    n_methods = len(all_method_keys)
    width = 0.22
    x = np.arange(n_methods)

    # Values for CellDiffusion and GCN
    values_cd = []
    values_gcn = []
    for method_key in all_method_keys:
        values_cd.append(scores_data.get((method_key, 'CellDiffusion'), {}).get(score_name, 0.0))
        values_gcn.append(scores_data.get((method_key, 'GCN'), {}).get(score_name, 0.0))

    cd_color, gcn_color = colors
    ax.bar(x - width / 2, values_cd, width,
        color=cd_color, edgecolor='black', linewidth=0.8, label='CellDiffusion')
    ax.bar(x + width / 2, values_gcn, width,
        color=gcn_color, edgecolor='black', linewidth=0.8, label='GCN')

    ax.set_xlabel('Graph Building Method', fontweight='bold')
    ax.set_ylabel('Score Value', fontweight='bold')
    ax.set_title(f'SCIB Evaluation: {score_name}', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(all_method_keys, rotation=45, ha='right')
    ax.legend(loc='best', framealpha=0.9)

    all_values = values_cd + values_gcn
    min_val = min(all_values) if all_values else 0
    max_val = max(all_values) if all_values else 1
    ax.set_ylim(bottom=max(0, min_val - 0.1), top=min(1.0, max_val + 0.1))

    # Add value labels
    for k, v in enumerate(values_cd):
        if v > 0.05:
            ax.text(k - width / 2, v + 0.01, f'{v:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    for k, v in enumerate(values_gcn):
        if v > 0.05:
            ax.text(k + width / 2, v + 0.01, f'{v:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    plt.tight_layout()

    print(f"Saving plot: {output_pdf}")
    Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
    plt.close()


# Plot and save three PDFs
score_colors = {
    'Total': ('#377eb8', '#a6cee3'),
    'Batch correction': ('#4daf4a', '#b2df8a'),
    'Bio conservation': ('#ff7f00', '#fdbf6f')
}

plot_single_score('Total', output_total, score_colors['Total'])
plot_single_score('Batch correction', output_batch, score_colors['Batch correction'])
plot_single_score('Bio conservation', output_bio, score_colors['Bio conservation'])

# Combined plot: Total and Bio conservation as subplots
print("\n=== Creating combined plot (Total + Bio conservation) ===")
fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)

def _plot_on_ax(ax, score_name: str, colors: tuple[str, str]) -> None:
    n_methods = len(all_method_keys)
    width = 0.22
    x = np.arange(n_methods)

    values_cd = [scores_data.get((m, 'CellDiffusion'), {}).get(score_name, 0.0) for m in all_method_keys]
    values_gcn = [scores_data.get((m, 'GCN'), {}).get(score_name, 0.0) for m in all_method_keys]

    cd_color, gcn_color = colors
    ax.bar(x - width / 2, values_cd, width,
        color=cd_color, edgecolor='black', linewidth=0.8, label='CellDiffusion')
    ax.bar(x + width / 2, values_gcn, width,
        color=gcn_color, edgecolor='black', linewidth=0.8, label='GCN')

    ax.set_title(score_name, fontsize=15, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(all_method_keys, rotation=45, ha='right', fontsize=11)

    all_values = values_cd + values_gcn
    min_val = min(all_values) if all_values else 0
    max_val = max(all_values) if all_values else 1
    ax.set_ylim(bottom=max(0, min_val - 0.1), top=min(1.0, max_val + 0.1))

    for k, v in enumerate(values_cd):
        if v > 0.05:
            ax.text(k - width / 2, v + 0.01, f'{v:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    for k, v in enumerate(values_gcn):
        if v > 0.05:
            ax.text(k + width / 2, v + 0.01, f'{v:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

_plot_on_ax(axes[0], 'Total', score_colors['Total'])
_plot_on_ax(axes[1], 'Bio conservation', score_colors['Bio conservation'])

axes[0].set_ylabel('Score Value', fontsize=14, fontweight='bold')
fig.suptitle('SCIB Evaluation: Total & Bio conservation', fontsize=16, fontweight='bold')

# Add separate legends for each subplot, positioned below each axis
axes[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), ncol=2, fontsize=11, framealpha=0.9)
axes[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), ncol=2, fontsize=11, framealpha=0.9)

plt.tight_layout(rect=[0, 0.25, 1, 0.95])

print(f"Saving combined plot: {output_combined}")
Path(output_combined).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_combined, dpi=300, bbox_inches='tight')
plt.close()

print("\n=== Plotting complete! ===")
print(f"  Graph methods: {all_method_keys}")
print(f"  Integration types: {integration_types}")
print(f"  Output PDFs: {output_total}, {output_batch}, {output_bio}, {output_combined}")
