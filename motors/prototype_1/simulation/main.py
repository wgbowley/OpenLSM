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

# Extract stripped values (remove units)
time = results.time.stripped
force = results.force.stripped
velocity = results.velocity.stripped
displacement = results.displacement.stripped
target_displacement = results.set_points.stripped
d_current = results.d_current.stripped
q_current = results.q_current.stripped
k_m = results.k_m.stripped
p_loss = results.p_loss.stripped
slot_temp = results.slot_temperature.stripped
slot_temp = results.slot_temperature.stripped
pole_temp = results.pole_temperature.stripped

plt.style.use("dark_background")
fig, axes = plt.subplots(5, 1, figsize=(12, 10), sharex=True)

fig.canvas.manager.set_window_title("Motor System Simulation Results")

MAIN_LW = 2.2
SEC_LW = 1.6
ALPHA = 0.85

# -------- Force & Currents --------
ax = axes[0]
ax.plot(time, force, color="#FF4500", linewidth=MAIN_LW, label="Force")
ax.set_ylabel("Force (N)", color="#FF4500")
ax.tick_params(colors="#FF4500")
ax.set_title("Force and Phase Currents vs Time", color="white")

ax_t = ax.twinx()
ax_t.plot(time, d_current, color="#00BFFF", linestyle="--",
            linewidth=SEC_LW, alpha=ALPHA, label="d-axis")
ax_t.plot(time, q_current, color="#7CFC00", linestyle="--",
            linewidth=SEC_LW, alpha=ALPHA, label="q-axis")
ax_t.set_ylabel("Current (A)", color="#00BFFF")
ax_t.tick_params(colors="#00BFFF")

ax.legend(loc="upper left", frameon=False)
ax_t.legend(loc="upper right", frameon=False)

# -------- Velocity & Displacement --------
ax = axes[1]
ax.plot(time, velocity, color="#00BFFF", linewidth=MAIN_LW, label="Velocity")
ax.set_ylabel("Velocity (m/s)", color="#00BFFF")
ax.tick_params(colors="#00BFFF")
ax.set_title("Velocity and Displacement vs Time", color="white")

ax_t = ax.twinx()
ax_t.plot(time, displacement, color="#FFA500",
            linewidth=MAIN_LW, label="Displacement")
ax_t.plot(time, target_displacement, color="#7CFC00",
            linestyle="--", linewidth=SEC_LW, label="Target")
ax_t.set_ylabel("Displacement (m)", color="#FFA500")
ax_t.tick_params(colors="#FFA500")

ax.legend(loc="upper left", frameon=False)
ax_t.legend(loc="upper right", frameon=False)

# -------- Motor Integrity --------
ax = axes[2]
ax.plot(time, k_m, color="#00BFFF", linewidth=MAIN_LW, label="Km")
ax.set_ylabel("Km (N/√W)", color="#00BFFF")
ax.tick_params(colors="#00BFFF")
ax.set_title("Motor Constant and Temperature", color="white")

ax_t = ax.twinx()
ax_t.plot(time, slot_temp, color="#FF4500",
            linewidth=SEC_LW, label="Slot Temp")
ax_t.set_ylabel("Temperature (K)", color="#FF4500")
ax_t.tick_params(colors="#FF4500")

ax.legend(loc="upper left", frameon=False)
ax_t.legend(loc="upper right", frameon=False)

# -------- Power Loss --------
ax = axes[3]
ax.plot(time, p_loss, color="#FF0000", linewidth=MAIN_LW)
ax.fill_between(time, p_loss, alpha=0.2)
ax.set_ylabel("Loss (W)", color="white")
ax.set_title("Resistive Power Loss", color="white")
ax.tick_params(colors="white")
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

# -------- Time Reference --------
axes[4].plot(time, time, color="#888888", alpha=0.0)  # dummy for spacing
axes[4].set_xlabel("Time (s)", color="white")
axes[4].set_title("Time Base", color="white")
axes[4].tick_params(colors="white")

fig.tight_layout(pad=2.0)
plt.show()