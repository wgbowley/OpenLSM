"""
Filename: main.py

Description:
    1D analytical simulation for
    a tubular linear synchronous motor. 
"""

from pathlib import Path
from picounits import Parser
from picounits import CURRENT, RESISTANCE, LENGTH, INDUCTANCE, POWER, NULLSET

from matplotlib import pyplot as plt
from model.solver import Solver


# Loads unit system, material library & parameters
ROOT_DIR = Path(__file__).resolve().parents[0]

# Materials & Parameter files
parameters_path = ROOT_DIR / "parameters.uiv"
parameters = Parser.open(parameters_path, ROOT_DIR / "../derived.ut")

# Magnetic Solver for tubular linear motor
solver = Solver(parameters)

# Calculates the sampling domain & step sizes
z_sample = parameters.numerics.displacement
step_size = parameters.numerics.samples.displacement_size

steps = int(z_sample / step_size)
offset = - z_sample / 2

# Begins position vs force sampling
print("-" * 50)
print(f"Steps: {steps}, Sample Range: ({offset:.3f}, {-offset:.3f})")
print("-" * 50)

# List to store
displacement = []
force_result = []

for index in range(0, steps):
    # Calculates the armature position
    z_output = index * step_size + offset
    z_solver = z_output - solver.armature_offset * LENGTH

    # Updates the currents & calculates force
    force = solver.compute_force(z_solver)

    force_result.append(force)
    displacement.append(z_output)

    # Prints position every interval
    if index % parameters.numerics.output.interval.stripped == 0:
        print(f"Position: {z_output} | Force: {force}")


# Prints a few different system parameters
print("-" * 50)
print(f"Slot Turns:         {solver.slot_turns * NULLSET}")
print(f"Line Resistance:    {solver.l2l_resistance * RESISTANCE}")
print(f"Line Inductance:    {solver.l2l_inductance * INDUCTANCE}")
print(f"Line Current:       {solver.phase_rms_current * CURRENT} (RMS)")
print(f"System Losses:      {solver.average_losses * POWER}  (1st order)")
print("-" * 50)


# Plotting position Vs force
plt.figure(figsize=(10, 5))
plt.plot(displacement, force_result, color="black")

# Graph styling details
plt.title(f"Analytical FOC Model: Force vs Position @ {solver.phase_rms_current * CURRENT} (RMS)")
plt.xlabel("Position (m)")
plt.ylabel("Force (N)")

# Vertical lines at stator boundaries
stator = parameters.stator.length.stripped
stator_half = stator / 2

plt.axvline(x=-stator_half, color='red', linestyle='--', linewidth=2, label='Stator End')
plt.axvline(x=stator_half, color='red', linestyle='--', linewidth=2)

plt.legend()
plt.grid(True, linestyle=':', alpha=0.6)
plt.tight_layout()

# Render plot to screen
plt.show()
