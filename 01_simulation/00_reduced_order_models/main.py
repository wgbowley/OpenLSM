"""
Filename: main.py

Description:
    Reduced order model for tubular 
    linear synchronous motor. 
    
    NOTE:
    This simulation model is still under development and 
    has not yet been validated or completed. 

    Do not assume the resulting design variables are suitable 
    for fabrication or real-world use without further analysis.
"""

from pathlib import Path
from math import sin, pi
import numpy as np
import matplotlib.pyplot as plt

from picounits import Parser
from src.model import TubularMotor
from src.solver import MagneticSolver

# Loads unit system, material library & parameters
ROOT_DIR = Path(__file__).resolve().parents[2]

# Materials & Parameter files
parameters_path = ROOT_DIR / "01_simulation/00_reduced_order_models/parameters.uiv"
materials_path = ROOT_DIR / "01_simulation/00_reduced_order_models/lib/materials.uiv"

materials = Parser.open(materials_path, ROOT_DIR / "metric.ut")
parameters = Parser.open(parameters_path, ROOT_DIR / "metric.ut")


# Load/Constructs motor model amd solver
model = TubularMotor(parameters, materials)
solver = MagneticSolver(model)

pha, phb, phc = model.armature.phase_a, model.armature.phase_b, model.armature.phase_c

# Set time to a fixed value (e.g., t=0)
time = 0
pha.current = model.raw_line_voltage / pha.resistance * sin(time)
phb.current = model.raw_line_voltage / pha.resistance * sin(time + 2 * pi/3)
phc.current = model.raw_line_voltage / pha.resistance * sin(time + 4 * pi /3)

rms_current = pha.current ** 2 + phb.current ** 2 + phc.current ** 2
rms_current = (rms_current / 3) ** 1/2

print(f"I_rms: {rms_current:.3f} A, copper_loss: {rms_current ** 2 * pha.resistance} W")

# Create position array (in meters)
position_array = np.linspace(-0.300, 0.300, 1000)
force_array = []

for position in position_array:
    force = solver.compute_force(position, 0.001)
    force_array.append(force)


fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(position_array * 1000, force_array, 'b-', linewidth=2)

# Auto-scale with margin
force_min = min(force_array)
force_max = max(force_array)
margin = abs(force_max - force_min) * 0.1
ax.set_ylim(force_min - margin, force_max + margin)

ax.set_xlabel('Position (mm)', fontsize=12)
ax.set_ylabel('Force (N)', fontsize=12)
ax.set_title(f'Force vs Position at t={time:.2f} s', fontsize=14)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()