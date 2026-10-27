"""
Filename: main.py

Description:
    1D analytical simulation for  
    a potted tubular annulus heat-sink.
"""

from pathlib import Path
from picounits import Parser

import matplotlib.pyplot as plt
from model.solver import Solver

# Loads unit system, material library & parameters
ROOT_DIR = Path(__file__).resolve().parents[0]

# Materials & Parameter files
parameters_path = ROOT_DIR / "parameters.uiv"
parameters = Parser.open(parameters_path, ROOT_DIR / "../derived.ut")

# Computes the performance grid across spacing / thickness.
solver = Solver(parameters)
axial, radial, dissipation = solver.compute_grid_performance()

# Maps the performance heat-map against input parameters
fig, ax = plt.subplots(figsize=(9, 6))
im = ax.imshow(
    dissipation,
    origin="lower",
    aspect="auto",
    extent=[axial[0], axial[-1], radial[0], radial[-1]],
    interpolation="bilinear",
    cmap='inferno'
)
cbar = fig.colorbar(im, ax=ax)

# Graph details
cbar.set_label("Heat dissipation (W)")
ax.set_xlabel("Axial Spacing (m)")
ax.set_ylabel("Radial Thickness (m)")
ax.set_title("Heat dissipation over solution space")

# Render plot to screen
fig.tight_layout()
plt.show()
