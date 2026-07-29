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
from math import sqrt

from picounits import Parser
from picounits import POWER, CURRENT, INDUCTANCE, RESISTANCE, LENGTH

from matplotlib import pyplot as plt

from src.physics import dynamic
from src.model import TubularMotor
from src.solver import MagneticSolver


# Loads unit system, material library & parameters
ROOT_DIR = Path(__file__).resolve().parents[0]

# Materials & Parameter files
parameters_path = ROOT_DIR / "parameters.uiv"
materials_path = ROOT_DIR / "lib/materials.uiv"

materials = Parser.open(materials_path, ROOT_DIR / "lib/metric.ut")
parameters = Parser.open(parameters_path, ROOT_DIR / "lib/metric.ut")

# Load/Constructs Motor Model amd Magnetic Solver
model = TubularMotor(parameters, materials)
solver = MagneticSolver(model)

# Calculates the displacement step size based on the armature (0-> 2pi basically)
steps = parameters.numerics.displacement_steps.stripped
step_size = model.raw_stator_length / steps

# Calculates the peak current based on the line voltage & phase resistance
# Assumes phase A, phase B and phase C have same resistance

# Breaks the phases to be addressable
pha, phb, phc = model.armature.phase_a, model.armature.phase_b, model.armature.phase_c
peak_current = model.raw_line_voltage / pha.resistance

# Calculates the rms current & copper losses for print out
rms_current = peak_current / sqrt(2)
losses = rms_current ** 2 * pha.resistance

print("-" * 10, " Static Parameters ",  "-" * 10)
print(f"Line & Peak: {peak_current*CURRENT}, RMS: {rms_current*CURRENT}")
print(f"Ind: {pha.inductance * INDUCTANCE}, Res: {pha.resistance * RESISTANCE}")
print(f"Losses: {losses * POWER:.3f}")
print("-" * 10, " Dynamic Parameters ",  "-" * 9)
print(f"Steps: {parameters.numerics.displacement_steps}, {step_size * LENGTH}")

half_length = model.raw_stator_length / 2
print(f"Start: {-half_length * LENGTH}, End: {half_length * LENGTH}")
print("-" * 41)

# Loop for FOC force output
displacement_data = []
force_data = []
for index in range(0, steps):
    # Calculates the electrical angle of the armature (Assumes stator, armature alignment)
    z_pos = step_size * index - half_length
    electrical_angle = dynamic.electrical_angle(z_pos, model.raw_pole_pitch)

    # Calculates the 3-phase currents from Q, B frame to Alpha, beta and finally A, B, C
    ab_ref = dynamic.inverse_park_transform(0, peak_current, electrical_angle)
    a, b, c = dynamic.inverse_clarke_transform(ab_ref)

    # Inputs the current for that position into the solver & solves for force at z_pos
    pha.current = a
    phb.current = b
    phc.current = c
    force = solver.compute_force(z_pos, parameters.numerics.integration_step_size.stripped)


    # Saves data from loop
    force_data.append(force)
    displacement_data.append(z_pos)


# Scales the data by 1000 from meters to millimeters
position_mm = [z * 1000.0 for z in displacement_data]

plt.figure(figsize=(10, 5))
plt.plot(position_mm, force_data, label="Reduced Order Model", color="blue", linewidth=2)

# Graph styling details
plt.title("Reduced Order Model: Force vs Position", fontsize=12)
plt.xlabel("Position (mm)", fontsize=10)
plt.ylabel("Force (N)", fontsize=10)

plt.grid(True, linestyle=':', alpha=0.6)
plt.tight_layout()

# Render plot to screen
plt.show()
