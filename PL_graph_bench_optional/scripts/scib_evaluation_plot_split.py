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

# Science journal figure style: 180mm double-column, 300+ DPI, sans-serif
# Ref: https://www.science.org/content/page/science-information-authors
MM_TO_IN = 1 / 25.4
COL_SINGLE_MM = 88
COL_DOUBLE_MM = 180
FIG_DPI = 600  # Science prefers 300-600 DPI
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.75,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": FIG_DPI,
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

    # Single-column width (88 mm) for Science
    w_in = COL_SINGLE_MM * MM_TO_IN
    fig, ax = plt.subplots(figsize=(w_in, w_in * 0.6))

    n_methods = len(all_method_keys)
    width = 0.28
    x = np.arange(n_methods)

    values_cd = [scores_data.get((m, 'CellDiffusion'), {}).get(score_name, 0.0) for m in all_method_keys]
    values_gcn = [scores_data.get((m, 'GCN'), {}).get(score_name, 0.0) for m in all_method_keys]

    cd_color, gcn_color = colors
    ax.bar(x - width / 2, values_cd, width, color=cd_color, edgecolor='#333333', linewidth=0.5, label='CellDiffusion')
    ax.bar(x + width / 2, values_gcn, width, color=gcn_color, edgecolor='#333333', linewidth=0.5, label='GCN')

    ax.set_xlabel('Graph building method')
    ax.set_ylabel('Score')
    ax.set_title(score_name)
    ax.set_xticks(x)
    ax.set_xticklabels(all_method_keys, rotation=45, ha='right')
    ax.legend(loc='upper right', frameon=True, framealpha=1, edgecolor='#cccccc')

    all_vals = values_cd + values_gcn
    y_min = max(0, (min(all_vals) if all_vals else 0) - 0.08)
    y_max = min(1.05, (max(all_vals) if all_vals else 1) + 0.05)
    ax.set_ylim(y_min, y_max)

    for k, v in enumerate(values_cd):
        if v > 0.05:
            ax.text(k - width / 2, v + 0.015, f'{v:.2f}', ha='center', va='bottom', fontsize=6)
    for k, v in enumerate(values_gcn):
        if v > 0.05:
            ax.text(k + width / 2, v + 0.015, f'{v:.2f}', ha='center', va='bottom', fontsize=6)

    plt.tight_layout()
    print(f"Saving plot: {output_pdf}")
    Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_pdf, dpi=FIG_DPI, bbox_inches='tight')
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

# Combined plot: Total and Bio conservation (Science double-column, 180 mm)
print("\n=== Creating combined plot (Total + Bio conservation) ===")
w_in = COL_DOUBLE_MM * MM_TO_IN
h_in = w_in * 0.38  # ~68 mm height
fig, axes_plot = plt.subplots(1, 2, figsize=(w_in, h_in))

def _plot_on_ax(ax, score_name: str, colors: tuple[str, str], panel_label: str) -> None:
    n_methods = len(all_method_keys)
    width = 0.28
    x = np.arange(n_methods)

    values_cd = [scores_data.get((m, 'CellDiffusion'), {}).get(score_name, 0.0) for m in all_method_keys]
    values_gcn = [scores_data.get((m, 'GCN'), {}).get(score_name, 0.0) for m in all_method_keys]

    cd_color, gcn_color = colors
    ax.bar(x - width / 2, values_cd, width, color=cd_color, edgecolor='#333333', linewidth=0.5, label='CellDiffusion')
    ax.bar(x + width / 2, values_gcn, width, color=gcn_color, edgecolor='#333333', linewidth=0.5, label='GCN')

    ax.set_title(score_name)
    ax.set_xlabel('Graph building method')
    ax.set_xticks(x)
    ax.set_xticklabels(all_method_keys, rotation=45, ha='right')

    all_vals = values_cd + values_gcn
    y_min = max(0, (min(all_vals) if all_vals else 0) - 0.08)
    y_max = min(1.05, (max(all_vals) if all_vals else 1) + 0.05)
    ax.set_ylim(y_min, y_max)

    for k, v in enumerate(values_cd):
        if v > 0.05:
            ax.text(k - width / 2, v + 0.015, f'{v:.2f}', ha='center', va='bottom', fontsize=6)
    for k, v in enumerate(values_gcn):
        if v > 0.05:
            ax.text(k + width / 2, v + 0.015, f'{v:.2f}', ha='center', va='bottom', fontsize=6)

    # Panel label (a, b) for Science
    ax.text(-0.12, 1.02, panel_label, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom')

_plot_on_ax(axes_plot[0], 'Total', score_colors['Total'], 'a')
_plot_on_ax(axes_plot[1], 'Bio conservation', score_colors['Bio conservation'], 'b')

axes_plot[0].set_ylabel('Score')
axes_plot[0].legend(loc='upper right', frameon=True, framealpha=1, edgecolor='#cccccc')
axes_plot[1].legend(loc='upper right', frameon=True, framealpha=1, edgecolor='#cccccc')

plt.tight_layout()
print(f"Saving combined plot: {output_combined}")
Path(output_combined).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_combined, dpi=FIG_DPI, bbox_inches='tight')
plt.close()

print("\n=== Plotting complete! ===")
print(f"  Graph methods: {all_method_keys}")
print(f"  Integration types: {integration_types}")
print(f"  Output PDFs: {output_total}, {output_batch}, {output_bio}, {output_combined}")
