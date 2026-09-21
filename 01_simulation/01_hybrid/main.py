"""
Filename: field_viewer.py

Description:
    This file allows for the FEMM solution
    to be viewed directly within python
    
    Uses: 
    pip install ifemm
"""

from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np

from ifemm import Parser


# Imports the parser and parses the .ans file
ROOT_DIR = Path(__file__).resolve().parents[0]
data = Parser.open(ROOT_DIR / 'resources/model.ans')

# Get the B field
length_unit = data.length_unit
x, y, bx, by = data.b_field()

# Calculate magnitude
b_magnitude = np.sqrt(bx**2 + by**2)

# Plot the solution
fig, ax = plt.subplots(figsize=(10, 8))

# Plot magnitude as background
contour = ax.contourf(x, y, b_magnitude, levels=50)
cbar = plt.colorbar(contour, ax=ax)
cbar.set_label('|B| (T)', fontsize=12)

ax.set_xlabel(f'x ({length_unit})', fontsize=12)
ax.set_ylabel(f'y ({length_unit})', fontsize=12)
ax.set_title('Magnetic Flux Density |B| with Field Lines', fontsize=14)
ax.axis('equal')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
