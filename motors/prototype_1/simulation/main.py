"""
Filename: main.py

Description:
    Prototype 1 linear Motor simulation script refer
    to configuration.uiv for implementation details
    and parameters.
    
    NOTE: 
    If this script fails to run. It most likely means 
    that the pyfea (v0.1.0) api has been deprecated
    
    NOTE:
    Dependencies: tabulate, matplotlib, pyfea, picounits
"""
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

from modules.initial_setup import pre_simulation_setup, initial_state
from modules.launch_analysis import Launch

from pyfea.models.tubular_linear_motor.main import TubularLinearMotor
from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver
from pyfea.solver.femm.domains.thermostatic.solver import FEMMThermostaticSolver

# Defines configuration file path and solver output folder path
path_lib = "motors/prototype_1/simulation/configuration.uiv"
solver_folder = "motors/prototype_1/simulation/outputs"

# Defines the tubular linear motor via configuration parameters
TubularMotor = TubularLinearMotor(path_lib)

# Initializes FEA solvers (magnetostatic and thermostatic)
Magnetic = FEMMMagnetostaticSolver(solver_folder)
Thermal = FEMMThermostaticSolver(solver_folder)

# Solves for initial magnetic parameters and initial thermal state
initial_conditions = pre_simulation_setup(TubularMotor, Magnetic)
initial_state(TubularMotor, Thermal)
print(initial_conditions)

# Initializes the launch analysis class
analysis = Launch(TubularMotor, Magnetic, Thermal, initial_conditions)
_, results = analysis.run()

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Motor Launch Analysis Results', fontsize=16)

# Extract stripped values (remove units)
time = results.time.stripped
force = results.force.stripped
velocity = results.velocity.stripped
displacement = results.displacement.stripped
target_displacement = results.set_points
d_current = results.d_current.stripped
q_current = results.q_current.stripped
k_m = results.k_m.stripped
p_loss = results.p_loss.stripped
slot_temp = results.slot_temperature.stripped
slot_temp = results.slot_temperature.stripped
pole_temp = results.pole_temperature.stripped

# Plot 1: Force and Currents
plt.style.use('seaborn-v0_8-whitegrid')

fig, axes = plt.subplots(2, 2, sharex=True, figsize=(13, 9))
fig.suptitle('Motor System Dynamics Overview', fontsize=16, fontweight='bold')

# ---------------- Plot 1: Force and Currents ----------------
ax1 = axes[0, 0]
ax1.plot(time, force, 'b-', label='Force', linewidth=2)
ax1.axhline(0, color='k', linewidth=0.8, alpha=0.2)
ax1.set_ylabel('Force (N)', color='b')
ax1.tick_params(axis='y', labelcolor='b')

ax1_twin = ax1.twinx()
ax1_twin.plot(time, d_current, 'r--', label='d-axis Current', alpha=0.7)
ax1_twin.plot(time, q_current, 'g--', label='q-axis Current', alpha=0.7)
ax1_twin.set_ylabel('Current (A)', color='r')
ax1_twin.tick_params(axis='y', labelcolor='r', alpha=0.7)
ax1_twin.spines['right'].set_alpha(0.5)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1_twin.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2,
           loc='upper center', bbox_to_anchor=(0.5, -0.18),
           ncol=3, frameon=False)

ax1.set_title('Force and Phase Currents')

# ---------------- Plot 2: Velocity and Displacement ----------------
ax2 = axes[0, 1]
ax2.plot(time, velocity, color='tab:blue', label='Velocity', linewidth=2)
ax2.axhline(0, color='k', linewidth=0.8, alpha=0.2)
ax2.set_ylabel('Velocity (m/s)', color='tab:blue')
ax2.tick_params(axis='y', labelcolor='tab:blue')

ax2_twin = ax2.twinx()
ax2_twin.plot(time, displacement, color='orange', label='Actual Displacement', linewidth=2)
ax2_twin.plot(time, target_displacement.stripped, 'g--', label='Target Displacement', linewidth=1.5, alpha=0.8)
ax2_twin.set_ylabel('Displacement (m)', color='orange')
ax2_twin.tick_params(axis='y', labelcolor='orange', alpha=0.7)
ax2_twin.spines['right'].set_alpha(0.5)

lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax2_twin.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2,
           loc='upper center', bbox_to_anchor=(0.5, -0.18),
           ncol=2, frameon=False)

ax2.set_title('Velocity and Displacement')

# ---------------- Plot 3: Motor Constant and Temperature ----------------
ax3 = axes[1, 0]
ax3.plot(time, k_m, color='tab:blue', label='Motor Constant (Km)', linewidth=2)
ax3.set_ylabel('Km (N/√W)', color='tab:blue')
ax3.tick_params(axis='y', labelcolor='tab:blue')

ax3_twin = ax3.twinx()
ax3_twin.plot(time, slot_temp, color='orange', label='Slot Temperature', linewidth=1.5)
ax3_twin.set_ylabel('Temperature (K)', color='orange')
ax3_twin.tick_params(axis='y', labelcolor='orange', alpha=0.7)
ax3_twin.spines['right'].set_alpha(0.5)

lines1, labels1 = ax3.get_legend_handles_labels()
lines2, labels2 = ax3_twin.get_legend_handles_labels()
ax3.legend(lines1 + lines2, labels1 + labels2,
           loc='upper center', bbox_to_anchor=(0.5, -0.18),
           ncol=2, frameon=False)

ax3.set_title('Motor Integrity')

# ---------------- Plot 4: Power Loss ----------------
ax4 = axes[1, 1]
ax4.fill_between(time, p_loss, alpha=0.15)
ax4.plot(time, p_loss, label='Resistive Loss', linewidth=2)
ax4.set_ylabel('Loss (W)')
ax4.set_xlabel('Time (s)')
ax4.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax4.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax4.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18),
           frameon=False)

ax4.set_title('Instantaneous Power Dissipation')

# Shared x label only on bottom row
ax3.set_xlabel('Time (s)')
ax4.set_xlabel('Time (s)')

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()