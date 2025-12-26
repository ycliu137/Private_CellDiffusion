"""
Plot SCIB evaluation results: aggregate scores vs k_add
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import re

# Get input and output from snakemake
input_table = snakemake.input.table
output_pdf = snakemake.output.pdf

print(f"\n=== Loading SCIB results ===")
print(f"Input file: {input_table}")

# Read results table
df = pd.read_csv(input_table, index_col=0)
print(f"Results table shape: {df.shape}")
print(f"Columns: {list(df.columns)}")
print(f"\nFirst few rows:")
print(df.head())

# Check if there's a "Metric Type" row (some scib_metrics versions include this)
# If so, we need to filter for aggregate scores
METRIC_TYPE_ROW = "Metric Type"
AGGREGATE_SCORE_TYPE = "Aggregate score"

has_metric_type_row = METRIC_TYPE_ROW in df.index
if has_metric_type_row:
    print(f"\nFound '{METRIC_TYPE_ROW}' row, filtering for aggregate scores...")
    # Get columns where metric type is "Aggregate score"
    aggregate_cols = df.columns[df.loc[METRIC_TYPE_ROW] == AGGREGATE_SCORE_TYPE].tolist()
    print(f"  Aggregate score columns (from Metric Type row): {aggregate_cols}")
    # Remove the metric type row from data
    df = df.drop(METRIC_TYPE_ROW)
    # Filter to only aggregate score columns (should be: Batch correction, Bio conservation, Total)
    if aggregate_cols:
        df = df[aggregate_cols]
        print(f"  Filtered to {len(aggregate_cols)} aggregate score columns")
    else:
        print(f"  Warning: No columns marked as '{AGGREGATE_SCORE_TYPE}' found, will search by column names")
        # Fallback: search for the three aggregate score columns by name
        expected_cols = ['Batch correction', 'Bio conservation', 'Total']
        aggregate_cols = [col for col in df.columns if col in expected_cols]
        if aggregate_cols:
            df = df[aggregate_cols]
            print(f"  Found aggregate score columns by name: {aggregate_cols}")
else:
    print(f"\nNo '{METRIC_TYPE_ROW}' row found, will search for aggregate score columns by name")
    # Try to find the three aggregate score columns directly
    expected_cols = ['Batch correction', 'Bio conservation', 'Total']
    aggregate_cols = [col for col in df.columns if col in expected_cols]
    if aggregate_cols:
        df = df[aggregate_cols]
        print(f"  Found aggregate score columns: {aggregate_cols}")

# Extract k_add values from method names (index)
# The index (Embedding column) corresponds to adata.obsm[x_dif_keys]
# Method names are the keys from adata.obsm, e.g., 'X_dif_kadd0', 'X_dif_kadd10'
print(f"\n=== Extracting k_add values from method names ===")
print(f"  Note: Method names (index/Embedding) correspond to adata.obsm[x_dif_keys]")
k_add_values = []
for method_name in df.index:
    # Try to find kadd in the method name
    # The method_name corresponds to adata.obsm[x_dif_keys], e.g., 'X_dif_kadd0', 'X_dif_kadd10'
    if 'kadd' in method_name.lower():
        try:
            # Extract the number after 'kadd' (case insensitive)
            match = re.search(r'kadd(\d+)', method_name, re.IGNORECASE)
            if match:
                k_add = int(match.group(1))
                k_add_values.append(k_add)
                print(f"  {method_name} -> k_add={k_add}")
            else:
                # Fallback: split on 'kadd' and take the first number
                parts = re.split(r'kadd', method_name, flags=re.IGNORECASE)
                if len(parts) > 1:
                    k_add_str = re.search(r'\d+', parts[1])
                    if k_add_str:
                        k_add = int(k_add_str.group())
                        k_add_values.append(k_add)
                        print(f"  {method_name} -> k_add={k_add}")
                    else:
                        print(f"Warning: Could not extract k_add from method name: {method_name}")
                        k_add_values.append(None)
                else:
                    print(f"Warning: Could not extract k_add from method name: {method_name}")
                    k_add_values.append(None)
        except (ValueError, AttributeError) as e:
            print(f"Warning: Could not extract k_add from method name: {method_name}, error: {e}")
            k_add_values.append(None)
    else:
        # If no 'kadd' in the name, skip this row (it's not one of our k_add variants)
        print(f"Skipping method name (no 'kadd'): {method_name}")
        k_add_values.append(None)

# Add k_add column to dataframe
df['k_add'] = k_add_values

# Filter out rows where k_add couldn't be extracted
df = df[df['k_add'].notna()].copy()
df = df.sort_values('k_add')

print(f"\n=== Extracting aggregate scores ===")

# At this point, df should already contain only the aggregate score columns
# The three aggregate scores are: "Batch correction", "Bio conservation", "Total"
# Map them to standardized names for plotting

found_scores = {}
expected_scores = {
    'Total': ['Total'],
    'Batch correction': ['Batch correction'],
    'Bio conservation': ['Bio conservation']
}

print(f"  Available columns after filtering: {list(df.columns)}")

# Direct mapping - check if the exact column names exist
for score_type, expected_names in expected_scores.items():
    for col in df.columns:
        if col in expected_names:
            found_scores[score_type] = col
            print(f"  Found '{score_type}' as column: {col}")
            break
        # Also check case-insensitive match
        col_lower = col.lower()
        score_type_lower = score_type.lower()
        if score_type_lower in col_lower or col_lower in score_type_lower:
            if score_type not in found_scores:  # Only add if not already found
                found_scores[score_type] = col
                print(f"  Found '{score_type}' as column: {col} (case-insensitive match)")

# If we still don't have all three, try fuzzy matching
if len(found_scores) < 3:
    print(f"  Warning: Only found {len(found_scores)} aggregate scores, trying fuzzy matching...")
    for col in df.columns:
        col_lower = col.lower()
        if 'total' in col_lower and 'Total' not in found_scores:
            found_scores['Total'] = col
            print(f"  Found 'Total' as column: {col}")
        elif ('batch' in col_lower or 'correction' in col_lower) and 'Batch correction' not in found_scores:
            found_scores['Batch correction'] = col
            print(f"  Found 'Batch correction' as column: {col}")
        elif ('bio' in col_lower or 'conservation' in col_lower) and 'Bio conservation' not in found_scores:
            found_scores['Bio conservation'] = col
            print(f"  Found 'Bio conservation' as column: {col}")

aggregate_scores = list(found_scores.keys())

if len(aggregate_scores) == 0:
    raise ValueError("Could not find any aggregate scores in the results table. Please check column names.")

print(f"\n=== Plotting aggregate scores vs k_add ===")

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))

# Plot each aggregate score
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, Orange, Green
markers = ['o', 's', '^']
linestyles = ['-', '--', '-.']

for i, score_type in enumerate(['Total', 'Batch correction', 'Bio conservation']):
    if score_type in found_scores:
        col_name = found_scores[score_type]
        values = df[col_name].values
        k_adds = df['k_add'].values
        
        # Sort by k_add to ensure proper line plotting
        sort_idx = np.argsort(k_adds)
        k_adds_sorted = k_adds[sort_idx]
        values_sorted = values[sort_idx]
        
        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]
        linestyle = linestyles[i % len(linestyles)]
        
        ax.plot(k_adds_sorted, values_sorted, 
                marker=marker, 
                linestyle=linestyle,
                color=color,
                linewidth=2,
                markersize=8,
                label=score_type)

# Customize plot
ax.set_xlabel('k_add', fontsize=12, fontweight='bold')
ax.set_ylabel('Aggregate Score', fontsize=12, fontweight='bold')
ax.set_title('SCIB Evaluation: Aggregate Scores vs k_add', fontsize=14, fontweight='bold')
ax.legend(loc='best', fontsize=10, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xlim(left=-5, right=max(df['k_add'].values) + 5)

# Improve layout
plt.tight_layout()

# Save figure
print(f"\n=== Saving plot ===")
print(f"Output file: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved successfully!")
print(f"  Plotted {len(aggregate_scores)} aggregate scores:")
for score_type in aggregate_scores:
    print(f"    - {score_type}")

print("\n=== Plotting complete! ===")

