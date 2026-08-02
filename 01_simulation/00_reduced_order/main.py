"""
Filename: main.py

Description:
    Reduced order model for tubular linear synchronous motor. 
"""

from pathlib import Path
from picounits import Parser
from matplotlib import pyplot as plt

from model.solver import Solver
from model.physics import field_oriented_control


# Loads unit system, material library & parameters
ROOT_DIR = Path(__file__).resolve().parents[0]

# Materials & Parameter files
parameters_path = ROOT_DIR / "parameters.uiv"
parameters = Parser.open(parameters_path, ROOT_DIR / "metric.ut")

# Magnetic Solver for tubular linear motor
solver = Solver(parameters)

time = parameters.numerics.time.stripped

# Calculates the sampling domain & step sizes
z_sample = parameters.numerics.sampling_size.stripped
step_size = parameters.numerics.displacement_step_size.stripped

steps = int(z_sample / step_size)
offset = - z_sample / 2

# List to store
displacement_data = []
force = []
for index in range(0, steps):
    # Calculates the electrical angle of the armature (Assumes stator, armature alignment)
    z_pos = index * step_size + offset
    electrical_angle = field_oriented_control.electrical_angle(z_pos, solver.dipole_axial_length)
    solver.update_currents(electrical_angle)

    # Calculates force over z distance
    force.append(solver.compute_force(z_pos))
    displacement_data.append(z_pos)


plt.figure(figsize=(10, 5))
plt.plot(displacement_data, force, label="Reduced Order Model", color="blue", linewidth=2)

# Graph styling details
plt.title("Reduced Order Model: Force vs Position", fontsize=12)
plt.xlabel("Position (m)", fontsize=10)
plt.ylabel("Force (N)", fontsize=10)

plt.grid(True, linestyle=':', alpha=0.6)
plt.tight_layout()

# Render plot to screen
plt.show()
