"""
Filename: thermal_optimize.py

Description:
    Prototype 1 linear Motor optimization script refer to design.uiv
    for initial motor parameters.
    
    NOTE:
    Dependencies: pyfea, picounits, numpy, matplotlib, BayesianOptimization
    
    If this script fails to run. It most likely means that
    the pyfea (v0.1.0) api has been deprecated.                          
"""

from operator import attrgetter
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from bayes_opt import BayesianOptimization

from module.initial_setup import static_evaluation
from model.tubular import TubularLinearMotor

from pyfea import watt, millimeter as mm
from pyfea.solver.femm.domains.thermostatic.solver import FEMMThermostaticSolver
from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver
from pyfea.solver.solver_outputs import SolverOutputs, ThermalOptions

# Parameters
power_losses = 40 * watt

# Defines configuration file and solver output path
BASE_DIR = Path(__file__).parent.parent.parent  
path_lib = BASE_DIR / "prototype_1/preliminary_design/parameters.uiv"
solver_folder = BASE_DIR / "prototype_1/preliminary_design/outputs"
csv_output = BASE_DIR / "prototype_1/preliminary_design/outputs/motor.csv"

# Defines the tubular linear motor via configuration parameters
TubularMotor = TubularLinearMotor(path_lib)

# Initializes FEA thermostatic & magnetostatic solvers
Thermal = FEMMThermostaticSolver(solver_folder)
Magnetic = FEMMMagnetostaticSolver(solver_folder)

# Why simulate both -> Because I am low-key lazy (like 10s at the start. its nothing!).
static_results = static_evaluation(TubularMotor, Thermal, Magnetic)

# Parameters
volumetric_heat = power_losses / static_results.slot_volume
initial_volume = static_results.heat_sink_volume
atmospheric_temp = TubularMotor.params.thermal.atmospheric_temperature

p_bounds = {
    "sink_radial": (0.2, 10),
    "fin_radial": (0.2, 10),
    "fin_axial": (0.2, 10),
    "sink_spacing": (0.2, 10),
}

# Counters
scores = []
iterations = []
def trial(sink_radial: float, fin_radial: float, fin_axial: float, sink_spacing: float) -> float:
    """ Tries a new heat sink design to improve thermals within the motor """
    sr, fr, fa, ss = sink_radial * mm, fin_radial * mm, fin_axial * mm, sink_spacing * mm
    try:
        TubularMotor.update_heat_sink_parameters(sr, fr, fa, ss)
        outputs = SolverOutputs()
        outputs.add_thermal(TubularMotor.SLOT_ID, ThermalOptions.AVERAGE_TEMPERATURE)
        outputs.add_thermal(TubularMotor.POLE_ID, ThermalOptions.AVERAGE_TEMPERATURE)
        outputs.add_thermal(TubularMotor.HEAT_SINK_ID, ThermalOptions.VOLUME)
        
        # Builds the motor
        domain = TubularMotor.construct_domain(Thermal)
        Thermal.setup(domain, "thermal_optimization")
        Thermal.update_heat_source(TubularMotor.armature_slots_material, volumetric_heat)
        
        # Solves and extract parameters
        thermal_results = Thermal.solve(outputs)
        
        volume = attrgetter(f"element_{TubularMotor.HEAT_SINK_ID.value}.volume")(thermal_results).value
        pole = attrgetter(f"element_{TubularMotor.POLE_ID.value}.average_temperature")(thermal_results).value
        slot = attrgetter(f"element_{TubularMotor.SLOT_ID.value}.average_temperature")(thermal_results).value
        
        temp = atmospheric_temp.value
        score = slot / temp + 2 * pole / temp #  + volume / initial_volume.value
        return -abs(score)
        
    except:
        # Fail-case if geometry is un-buildable
        return -10
    
def objective(sink_radial, fin_radial, fin_axial, sink_spacing):
    """ Wrapper for main trial function """
    score = trial(sink_radial, fin_radial, fin_axial, sink_spacing)

    scores.append(score)
    iterations.append(len(scores))

    return score


optimizer = BayesianOptimization(
    f=objective,
    pbounds=p_bounds,
    verbose=2,
    random_state=1
)
optimizer.maximize(init_points=100, n_iter=200)

plt.figure()

plt.plot(iterations, scores)
plt.xlabel("Iteration")
plt.ylabel("Thermal Score")
plt.title("Thermal Optimization Progress")

plt.show()