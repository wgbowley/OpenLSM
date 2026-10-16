"""
Filename: search_tradeoffs.py

Description:
    Loads the csv trace file and plots the 
    different metrics against each other.
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file from path to memory
ROOT_DIR = Path(__file__).resolve().parents[0]
csv_filename = ROOT_DIR / "09_run_18_slots/trace.csv"

df = pd.read_csv(csv_filename)

metrics = [
    "Absolute Force (N)",
    "Force per Ampere (N/A)",
    "Force slew rate (N/s)",
    "Copper Loss (W)",
    "Wire diameter (m)",
    "Radial Thickness (m)",
]

n = len(metrics)
fig, axes = plt.subplots(n, n, figsize=(3.2 * n, 3.2 * n))

iters = df["Iteration"].to_numpy()

for r in range(n):
    # Plots the scatter plot for each metric vs each metric.
    for c in range(n):
        ax = axes[r, c]
        y_col = metrics[r]
        x_col = metrics[c]

        # Scatter
        ax.scatter(df[x_col], df[y_col])

        ax.set_xlabel(x_col, fontsize=7)
        ax.set_ylabel(y_col, fontsize=7)
        ax.tick_params(labelsize=7)


fig.tight_layout(pad=1.5, h_pad=1.2, w_pad=1.2)

out_path = csv_filename.parent / "tradeoffs.png"
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Visualization saved successfully as '{out_path}'.")

plt.show()
