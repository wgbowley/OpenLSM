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

from modules.initial_setup import pre_simulation_setup, initial_state
from modules.launch_analysis import Launch

from pyfea import dimensionless
from pyfea.models.tubular_linear_motor.main import TubularLinearMotor
from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver
from pyfea.solver.femm.domains.thermostatic.solver import FEMMThermostaticSolver

from pyfea.solver.solver_outputs import SolverOutputs, ThermalOptions


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
d_current = results.d_current.stripped
q_current = results.q_current.stripped
mechanical_energy = results.mechanical_energy.stripped
electrical_energy = results.electrical_energy.stripped
slot_temp = results.slot_temperature.stripped
pole_temp = results.pole_temperature.stripped

# Plot 1: Force and Currents
ax1 = axes[0, 0]
ax1.plot(time, force, 'b-', label='Force', linewidth=2)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Force (N)', color='b')
ax1.tick_params(axis='y', labelcolor='b')
ax1.grid(True, alpha=0.3)

ax1_twin = ax1.twinx()
ax1_twin.plot(time, d_current, 'r--', label='d-axis Current', alpha=0.7)
ax1_twin.plot(time, q_current, 'g--', label='q-axis Current', alpha=0.7)
ax1_twin.set_ylabel('Current (A)', color='r')
ax1_twin.tick_params(axis='y', labelcolor='r')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1_twin.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
ax1.set_title('Force and Phase Currents vs Time')

# Plot 2: Velocity and Displacement
ax2 = axes[0, 1]
ax2.plot(time, velocity, 'purple', label='Velocity', linewidth=2)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Velocity (m/s)', color='purple')
ax2.tick_params(axis='y', labelcolor='purple')
ax2.grid(True, alpha=0.3)

ax2_twin = ax2.twinx()
ax2_twin.plot(time, displacement, 'orange', label='Displacement', linewidth=2)
ax2_twin.set_ylabel('Displacement (m)', color='orange')
ax2_twin.tick_params(axis='y', labelcolor='orange')

lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax2_twin.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
ax2.set_title('Velocity and Displacement vs Time')

# Plot 3: Energy
ax3 = axes[1, 0]
ax3.plot(time, mechanical_energy, 'b-', label='Mechanical Energy', linewidth=2)
ax3.plot(time, electrical_energy, 'r-', label='Electrical Energy', linewidth=2)
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Energy (J)')
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_title('Energy vs Time')

# Plot 4: Temperatures
ax4 = axes[1, 1]
ax4.plot(time, slot_temp, 'r-', label='Slot Temperature', linewidth=2)
ax4.plot(time, pole_temp, 'b-', label='Pole Temperature', linewidth=2)
ax4.set_xlabel('Time (s)')
ax4.set_ylabel('Temperature (K)')
ax4.legend()
ax4.grid(True, alpha=0.3)
ax4.set_title('Temperature vs Time')

plt.tight_layout()
plt.show()