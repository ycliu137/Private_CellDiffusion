"""
Plot metrics from metrics_log.csv
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# Load input and output paths
input_csv = snakemake.input.metrics
output_pdf = snakemake.output.pdf

print(f"Reading metrics from: {input_csv}")
df = pd.read_csv(input_csv)

# Check required columns
required_cols = ["k_added", "neighbor_purity_before", "neighbor_purity_after"]
missing_cols = [col for col in required_cols if col not in df.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

# Sort by k_added to ensure proper line plotting
df = df.sort_values("k_added")

# Create the plot
print("Creating plot...")
fig, ax = plt.subplots(figsize=(8, 6))

# Plot both lines
ax.plot(df["k_added"], df["neighbor_purity_before"], 
        marker='o', label='Neighbor Purity Before', linewidth=2)
ax.plot(df["k_added"], df["neighbor_purity_after"], 
        marker='s', label='Neighbor Purity After', linewidth=2)

# Customize the plot
ax.set_xlabel('k_added', fontsize=12)
ax.set_ylabel('Neighbor Purity', fontsize=12)
ax.set_title('Neighbor Purity vs k_added', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(left=min(df["k_added"]) - 5, right=max(df["k_added"]) + 5)

# Tight layout for better appearance
plt.tight_layout()

# Save as PDF
print(f"Saving plot to: {output_pdf}")
Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
plt.close()

print("Plot saved successfully!")

