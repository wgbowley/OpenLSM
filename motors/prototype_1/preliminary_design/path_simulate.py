"""
Filename: simulate.py

Description:
    Prototype 1 linear Motor simulation script refer to parameters.uiv 
    for implementation details & parameters.
    
    NOTE: 
    Dependencies: matplotlib, pyfea, picounits
    
    If this script fails to run. It most likely means 
    that the pyfea (v0.1.0) api has been deprecated
"""

import matplotlib.pyplot as plt

from pathlib import Path
from module.sim_definitions import PathSegment
from module.initial_setup import static_evaluation
from module.dynamic_analysis import PointToPoint
from model.tubular import TubularLinearMotor

from pyfea import second, millimeter as mm, dimensionless as D
from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver
from pyfea.solver.femm.domains.thermostatic.solver import FEMMThermostaticSolver

# Defines configuration file and solver output path
BASE_DIR = Path(__file__).parent.parent.parent  
path_lib = BASE_DIR / "prototype_1/preliminary_design/design.uiv"
solver_folder = BASE_DIR / "prototype_1/preliminary_design/outputs"
csv_output = BASE_DIR / "prototype_1/preliminary_design/outputs/motor.csv"

# Defines the tubular linear motor via configuration parameters
TubularMotor = TubularLinearMotor(path_lib)

# Initializes FEA solvers (magnetostatic and thermostatic)
Magnetic = FEMMMagnetostaticSolver(solver_folder)
Thermal = FEMMThermostaticSolver(solver_folder)

static_results = static_evaluation(TubularMotor, Thermal, Magnetic)
print(static_results)

analysis = PointToPoint(TubularMotor, Thermal, Magnetic, static_results)
segment = PathSegment(50 * mm, 200 * mm/second, 10000 * mm/second**2, 0.5 * second)

results, summary = analysis.run(segment, csv_output, True)
print(summary)

# Extract stripped values (removes units)
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

# Setup for 3 clean graphs
fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
fig.suptitle("MOTOR SYSTEM PERFORMANCE SIMULATION", fontsize=16, fontweight='bold', color='white')

# Style Constants
MAIN_LW = 2.0
SEC_LW = 1.5
ALPHA = 0.3

# 1. DYNAMICS: Force & Phase Currents
ax1 = axes[0]
ax1.plot(time, force, color="#00FFD1", linewidth=MAIN_LW, label="Output Force")
ax1.set_ylabel("FORCE (N)", fontweight='bold')

ax1_t = ax1.twinx()
ax1_t.plot(time, d_current, color="#FF00FF", linestyle="--", linewidth=SEC_LW, label="$I_d$ (d-axis)")
ax1_t.plot(time, q_current, color="#FFFF00", linestyle="--", linewidth=SEC_LW, label="$I_q$ (q-axis)")
ax1_t.set_ylabel("CURRENTS (A)", fontweight='bold')

ax1.set_title("Electromechanical Dynamics", loc='left', alpha=0.7)
ax1.legend(loc="upper left", frameon=False)
ax1_t.legend(loc="upper right", frameon=False)

# 2. KINEMATICS: Velocity & Displacement
ax2 = axes[1]
ax2.plot(time, velocity, color="#00BFFF", linewidth=MAIN_LW, label="Velocity")
ax2.set_ylabel("VELOCITY (m/s)", fontweight='bold')

ax2_t = ax2.twinx()
ax2_t.plot(time, displacement, color="#FF4500", linewidth=MAIN_LW, label="Actual Pos")
ax2_t.plot(time, target_displacement, color="#7CFC00", linestyle=":", linewidth=SEC_LW, label="Target")
ax2_t.set_ylabel("DISPLACEMENT (m)", fontweight='bold')

ax2.set_title("Motion Tracking", loc='left', alpha=0.7)
ax2.legend(loc="upper left", frameon=False)
ax2_t.legend(loc="upper right", frameon=False)

# 3. THERMAL & LOSS: Motor Integrity & Power
ax3 = axes[2]
ax3.plot(time, slot_temp, color="#FF3131", linewidth=MAIN_LW, label="Slot Temp")
ax3.fill_between(time, slot_temp, alpha=0.1, color="#FF3131")
ax3.set_ylabel("TEMPERATURE (K)", fontweight='bold')

ax3_t = ax3.twinx()
ax3_t.plot(time, p_loss, color="#FFFFFF", linewidth=SEC_LW, alpha=0.6, label="Power Loss")
ax3_t.set_ylabel("LOSS (W)", fontweight='bold')

ax3.set_title("Thermal Load & Efficiency", loc='left', alpha=0.7)
ax3.set_xlabel("TIME", fontweight='bold')
ax3.legend(loc="upper left", frameon=False)
ax3_t.legend(loc="upper right", frameon=False)

for i, ax in enumerate(axes):
    # 1. Enable Y-axis ticks and labels for the main axes
    ax.tick_params(axis='y', colors='white', labelsize=9)
    ax.yaxis.set_visible(True)
    
    # 2. Find and enable labels for the twin axes (the ones on the right)
    # We look through all axes in the figure that aren't the main ones
    for other_ax in fig.axes:
        if other_ax not in axes:
            other_ax.tick_params(axis='y', colors='white', labelsize=9)
            other_ax.yaxis.set_visible(True)

    # 3. Enable X-axis labels only for the bottom plot (since sharex=True)
    if i == 2:
        ax.tick_params(axis='x', colors='white', labelsize=9)
        ax.xaxis.set_visible(True)

    # Aesthetic styling (Keeping the clean look but showing the lines)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#444444')
    ax.spines['bottom'].set_color('#444444')
    ax.grid(True, linestyle=':', alpha=0.2)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()