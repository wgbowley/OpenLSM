"""
Filename: main.py

Description:
    Hybrid simulation for a tubular 
    linear synchronous motor.
"""

from pathlib import Path
import numpy as np

from ifemm import Parser
from picounits import ENERGY


# Imports the parser and parses the .ans file
ROOT_DIR = Path(__file__).resolve().parents[0]
data = Parser.open(ROOT_DIR / 'resources/model.ans')

# Get the B field
length_unit = data.length_unit
scale = data.length_scale
x, y, bx, by = data.b_field()

# Calculate b-field magnitude & removes invalid entries.
b_magnitude = np.sqrt(bx**2 + by**2)
valid = np.isfinite(b_magnitude)

mu = 4 * np.pi * 1e-7

# Ensures both the x & y grids are arrays
x = np.asarray(x)
y = np.asarray(y)

# Assumes uniform x & y grid
dx = np.abs(x[0, 1] - x[0, 0]) * scale
dy = np.abs(y[1, 0] - y[0, 0]) * scale
r = x * scale

# Integrate only the valid cells
# Axisymmetric volume element: dV = 2*pi*r*dr*dz
cell_volume = 2 * np.pi * np.abs(r) * dx * dy
U = np.sum(b_magnitude[valid]**2 * cell_volume[valid]) / (2 * mu)

print(f"Stored magnetic energy U = {U*ENERGY:.3f}")
