"""
Plot neighbor purity and bio conservation score in the same bar plot
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np

# Get input and output from snakemake
input_purity_celldiffusion = snakemake.input.purity_celldiffusion
input_purity_gcn = snakemake.input.purity_gcn
input_scib_table = snakemake.input.scib_table
output_pdf = snakemake.output.pdf

print(f"\n=== Loading data ===")
print(f"Purity CSV (CellDiffusion): {input_purity_celldiffusion}")
print(f"Purity CSV (GCN): {input_purity_gcn}")
print(f"SCIB results table: {input_scib_table}")

# Load purity results
purity_cd = pd.read_csv(input_purity_celldiffusion)
purity_gcn = pd.read_csv(input_purity_gcn)

print(f"\n=== Purity results ===")
print(f"CellDiffusion purity: {purity_cd['neighbor_purity'].iloc[0]:.4f}")
print(f"GCN purity: {purity_gcn['neighbor_purity'].iloc[0]:.4f}")

# Load SCIB results table
scib_df = pd.read_csv(input_scib_table, index_col=0)
print(f"\nSCIB results table shape: {scib_df.shape}")
print(f"Columns: {list(scib_df.columns)}")
print(f"Index (Embedding methods): {list(scib_df.index)}")

# Check for "Metric Type" row and remove it if present
METRIC_TYPE_ROW = "Metric Type"
if METRIC_TYPE_ROW in scib_df.index:
    print(f"\nRemoving '{METRIC_TYPE_ROW}' row from data")
    scib_df = scib_df.drop(METRIC_TYPE_ROW)

# Find CellDiffusion (X_dif) and GCN (X_gcn) rows
celldiffusion_row = None
gcn_row = None

print(f"\n=== Finding methods in SCIB results ===")
for method_name in scib_df.index:
    method_str = str(method_name)
    if 'X_dif' in method_str and 'X_dif_kadd' not in method_str:
        celldiffusion_row = method_name
        print(f"Found CellDiffusion: {method_name}")
    elif 'X_gcn' in method_str:
        gcn_row = method_name
        print(f"Found GCN: {method_name}")

if celldiffusion_row is None:
    raise ValueError("Could not find CellDiffusion (X_dif) in SCIB results table. Available methods: " + str(list(scib_df.index)))

if gcn_row is None:
    raise ValueError("Could not find GCN (X_gcn) in SCIB results table. Available methods: " + str(list(scib_df.index)))

# Extract Bio conservation score
bio_conservation_col = None
print(f"\n=== Extracting Bio conservation score ===")
for col in scib_df.columns:
    col_lower = col.lower()
    if ('bio' in col_lower and 'conservation' in col_lower) or ('conservation' in col_lower and 'bio' in col_lower):
        bio_conservation_col = col
        print(f"Found Bio conservation column: {col}")
        break

if bio_conservation_col is None:
    raise ValueError("Could not find Bio conservation column. Available columns: " + str(list(scib_df.columns)))

# Convert to numeric
scib_df[bio_conservation_col] = pd.to_numeric(scib_df[bio_conservation_col], errors='coerce')

# Extract bio conservation scores
bio_conservation_cd = float(scib_df.loc[celldiffusion_row, bio_conservation_col])
bio_conservation_gcn = float(scib_df.loc[gcn_row, bio_conservation_col])

if pd.isna(bio_conservation_cd):
    print(f"  Warning: Bio conservation is NaN for CellDiffusion, setting to 0.0")
    bio_conservation_cd = 0.0
if pd.isna(bio_conservation_gcn):
    print(f"  Warning: Bio conservation is NaN for GCN, setting to 0.0")
    bio_conservation_gcn = 0.0

print(f"CellDiffusion Bio conservation: {bio_conservation_cd:.4f}")
print(f"GCN Bio conservation: {bio_conservation_gcn:.4f}")

# Prepare data for plotting
purity_cd_val = purity_cd['neighbor_purity'].iloc[0]
purity_gcn_val = purity_gcn['neighbor_purity'].iloc[0]

methods = ['CellDiffusion', 'GCN']
metrics = ['Neighbor Purity', 'Bio Conservation']
data = {
    'CellDiffusion': [purity_cd_val, bio_conservation_cd],
    'GCN': [purity_gcn_val, bio_conservation_gcn]
}

# Create bar plot
print(f"\n=== Creating bar plot ===")

fig, ax = plt.subplots(figsize=(10, 6))

# Set up bar positions
x = np.arange(len(metrics))
width = 0.35  # Width of bars

colors = ['#1f77b4', '#ff7f0e']  # Blue for CellDiffusion, Orange for GCN

# Plot bars
for i, method in enumerate(methods):
    values = data[method]
    offset = (i - 0.5) * width
    ax.bar(x + offset, values, width, label=method, color=colors[i], alpha=0.8, edgecolor='black', linewidth=1.2)

# Customize plot
ax.set_xlabel('Metric', fontsize=12, fontweight='bold')
ax.set_ylabel('Score Value', fontsize=12, fontweight='bold')
ax.set_title('Neighbor Purity and Bio Conservation: CellDiffusion vs GCN', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(metrics, rotation=0, ha='center')
ax.legend(loc='best', fontsize=11, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--', axis='y')

# Set y-axis limits dynamically based on max value
all_values = []
for method in methods:
    all_values.extend(data[method])
max_value = max(all_values) if all_values else 1.0
y_max = min(1.0, max_value * 1.1)  # Add 10% padding, cap at 1.0
ax.set_ylim(bottom=0, top=y_max)

# Add value labels on bars
for i, method in enumerate(methods):
    values = data[method]
    offset = (i - 0.5) * width
    for j, v in enumerate(values):
        ax.text(j + offset, v + 0.01, f'{v:.3f}', 
                ha='center', va='bottom', fontsize=9, fontweight='bold')

plt.tight_layout()

# Save figure
print(f"\n=== Saving plot ===")
print(f"Output file: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print(f"Bar plot saved successfully!")
print(f"  Compared methods: {methods}")
print(f"  Metrics: {metrics}")

print("\n=== Plotting complete! ===")

