"""
Filename: visualize.py

Description:
    Visualizes the normalized score improvement with 
    respect to slot thickness and slot axial length.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

# Load the CSV file from path to memory
ROOT_DIR = Path(__file__).resolve().parents[0]
csv_filename = ROOT_DIR / "08_run_300um_18_slots/trace.csv"
df = pd.read_csv(csv_filename)


best = []
current_max = -float("inf")
for idx, row in df.iterrows():
    # Filters for best normalized score over iterations.
    score = row["Normalized Score"]

    if score > current_max:
        current_max = score
        best.append(idx)

df_best = df.loc[best].copy()

x = df_best["Axial Length (m)"]
y = df_best["Radial Thickness (m)"]
scores = df_best["Normalized Score"]
iterations = df_best["Iteration"]

# The absolute best is the last row in df_best
best_row = df_best.iloc[-1]
best_iter = int(best_row["Iteration"])

plt.figure(figsize=(9, 7))
ax = plt.gca()

# Scatters all of the points
x_all = df["Axial Length (m)"]
y_all = df["Radial Thickness (m)"]
scores_all = df["Normalized Score"]

scatter_all = ax.scatter(x_all, y_all, c=scores_all)

for _, row in df_best.iterrows():
    bx = row["Axial Length (m)"]
    by = row["Radial Thickness (m)"]
    step = row["step_size (m)"]

    # Draws a box around the current search area.
    rect = patches.Rectangle(
        (bx - step, by - step),
        2 * step,
        2 * step,
        linewidth=1.5,
        edgecolor="black",
        facecolor="none",
        linestyle="--",
        alpha=0.6,
        zorder=2,
    )
    ax.add_patch(rect)


# Draws trajectory line connecting best at-time results
plt.plot(x, y, color="tab:blue", linestyle="-", linewidth=1, label="Path")
plt.xlim(xmin=0)
plt.ylim(ymin=0)

cbar = plt.colorbar(scatter_all)
cbar.set_label("Normalized Score", fontsize=12)

ax.set_aspect("equal", adjustable="datalim")

# Labels and styling
plt.xlabel("Axial Length (m)")
plt.ylabel("Radial Thickness (m)")
plt.title("Stochastic local search with adaptive step-size reduction")
plt.grid(True, linestyle="--")
plt.legend(loc="best")

plt.tight_layout()
plt.savefig("trace.png", dpi=300)
print("Visualization saved successfully as 'trace.png'.")
plt.show()
