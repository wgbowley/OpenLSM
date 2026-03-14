"""
Filename: static_simulate.py

Description:
    Prototype 1 linear Motor simulation script refer to design.uiv
    for motor parameters.
    
    NOTE:
    Dependencies: pyfea
    
    If this script fails to run. It most likely means that
    the pyfea (v0.1.0) api has been deprecated.                          
"""


from operator import attrgetter
from pathlib import Path

from pyfea import watt
from model.tubular import TubularLinearMotor
from module.initial_setup import static_evaluation

from pyfea.solver.femm.domains.thermostatic.solver import FEMMThermostaticSolver
from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver
from pyfea.solver.solver_outputs import SolverOutputs, ThermalOptions

# Defines configuration file and solver output path
BASE_DIR = Path(__file__).parent.parent.parent  
path_lib = BASE_DIR / "prototype_1/preliminary_design/parameters.uiv"
solver_folder = BASE_DIR / "prototype_1/preliminary_design/outputs"
csv_output = BASE_DIR / "prototype_1/preliminary_design/outputs/motor.csv"

# Heat source
power_losses = 54 * watt

# Defines the tubular linear motor via configuration parameters
TubularMotor = TubularLinearMotor(path_lib)
Magnetic = FEMMMagnetostaticSolver(solver_folder)

# Initializes FEA thermostatic & magnetostatic solvers
Thermal = FEMMThermostaticSolver(solver_folder)

# Why simulate both -> Because I am low-key lazy (like 10s at the start. its nothing!).
static_results = static_evaluation(TubularMotor, Thermal, Magnetic)
print(static_results)

# Parameters
volumetric_heat = power_losses / static_results.slot_volume
initial_volume = static_results.heat_sink_volume

outputs = SolverOutputs()
outputs.add_thermal(TubularMotor.SLOT_ID, ThermalOptions.AVERAGE_TEMPERATURE)
outputs.add_thermal(TubularMotor.POLE_ID, ThermalOptions.AVERAGE_TEMPERATURE)
outputs.add_thermal(TubularMotor.HEAT_SINK_ID, ThermalOptions.VOLUME)

# Builds the motor
domain = TubularMotor.construct_domain(Thermal)
Thermal.setup(domain)
Thermal.update_heat_source(TubularMotor.armature_slots_material, volumetric_heat)

# Solves and extract parameters
thermal_results = Thermal.solve(outputs)

volume = attrgetter(f"element_{TubularMotor.HEAT_SINK_ID.value}.volume")(thermal_results)
pole = attrgetter(f"element_{TubularMotor.POLE_ID.value}.average_temperature")(thermal_results)
slot = attrgetter(f"element_{TubularMotor.SLOT_ID.value}.average_temperature")(thermal_results)

cost = TubularMotor.params.material_cost.aluminum
mass = volume * TubularMotor.heat_sink_material.values().physical.density

print(f"loss: {power_losses:.3f}, pole_temp: {pole:.3f}, slot_temp: {slot:.3f}, heat_sink: ${mass*cost:.3f}")