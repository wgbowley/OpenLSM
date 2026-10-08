"""
Filename: main.py

Description:
    Hybrid simulation for a tubular 
    linear synchronous motor.
"""

from math import sin, pi
from pathlib import Path
from matplotlib import pyplot as plt

from ifemm import Parser as iParser
from picounits import Parser

from model.solver import Solver

# Imports the parser and parses the .ans file
ROOT_DIR = Path(__file__).resolve().parents[0]
data = iParser.open(ROOT_DIR / 'resources/model.ans')

# Materials & Parameter files
parameters_path = ROOT_DIR / "parameters.uiv"
parameters = Parser.open(parameters_path, ROOT_DIR / "../derived.ut")

solver = Solver(parameters, data)

solver.i_pha, solver.i_phb, solver.i_phc = 10 * sin(0), 10 * sin(2 * pi / 3), 10 * sin(4 * pi / 3)


r = solver.r_eval          # (n_r, n_z, 3)
B = solver.b_field         # (n_r, n_z)

R = r[..., 0]              # (n_r, n_z) — x-coordinate = radial
Z = r[..., 2]              # (n_r, n_z) — z-coordinate

fig, ax = plt.subplots(figsize=(8, 5))
im = ax.pcolormesh(R, Z, B, shading="auto", cmap="viridis")
ax.set_xlabel("r")
ax.set_ylabel("z")
ax.set_aspect("equal")
plt.colorbar(im, label="|B|")
plt.title("Kernel field magnitude")
plt.show()