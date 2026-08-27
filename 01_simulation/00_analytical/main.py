"""
Filename: main.py

Description:
    1D analytical simulation for
    a tubular linear synchronous motor. 
"""

from pathlib import Path
from picounits import CURRENT, RESISTANCE, LENGTH, FORCE, INDUCTANCE, POWER, NULLSET, Parser
from matplotlib import pyplot as plt

from model.solver import Solver


# Loads unit system, material library & parameters
ROOT_DIR = Path(__file__).resolve().parents[0]

# Materials & Parameter files
parameters_path = ROOT_DIR / "parameters.uiv"
parameters = Parser.open(parameters_path, ROOT_DIR / "../metric.ut")

# Magnetic Solver for tubular linear motor
solver = Solver(parameters)


# Calculates the sampling domain & step sizes
z_sample = parameters.numerics.displacement.stripped
step_size = parameters.numerics.samples.displacement_size.stripped

steps = int(z_sample / step_size)
offset = - z_sample / 2

# Begins position vs force sampling
print(f"Steps: {steps}, Sample Range: ({offset * LENGTH:.3f}, {-offset * LENGTH:.3f})")
print("-" * 20)

# List to store
displacement = []
force_magnitude = []

for index in range(0, steps):
    # Calculates the armature position
    z_output = index * step_size + offset
    z_solver = z_output - solver.armature_offset

    # Updates the currents & calculates force
    force = solver.compute_force(z_solver)

    # Calculates force over z distance
    force_magnitude.append(force)
    displacement.append(z_output)

    # Prints position every interval
    if index % parameters.numerics.output.interval.stripped == 0:
        print(f"Position: {z_output * LENGTH:.3f} | Force: {force * FORCE:.3f}")


# Prints a few different system parameters
print("-" * 20)
print(f"Slot Turns:         {solver.slot_turns * NULLSET:.3f}")
print(f"Line Resistance:    {solver.l2l_resistance * RESISTANCE:.3f}")
print(f"Line Inductance:    {solver.l2l_inductance * INDUCTANCE:.3f}")
print(f"Line Current:       {solver.phase_rms_current * CURRENT:.3f} (RMS)")
print(f"System Losses:      {solver.average_losses * POWER:.3f} (Copper)")
print("-" * 20)


# Plotting position Vs magnitude of force
plt.figure(figsize=(10, 5))
plt.plot(displacement, force_magnitude, color="black")

# Graph styling details
plt.title("Analytical Model: Force vs Position")
plt.xlabel("Position (m)")
plt.ylabel("|Force| (N)")

plt.grid(True, linestyle=':', alpha=0.6)
plt.tight_layout()

# Render plot to screen
plt.show()
