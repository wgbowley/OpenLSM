"""
Filename: main.py

Description:
    Prototype 1 linear Motor simulation / optimization 
    script refer to configuration.uiv for implementation 
    details and parameters.
    
    NOTE: 
    Dependencies: matplotlib, pyfea, picounits, pymoo
    
    If this script fails to run. It most likely means 
    that the pyfea (v0.1.0) api has been deprecated
    
"""

from module.initial_setup import static_evaluation

from pyfea.models.tubular_linear_motor.main import TubularLinearMotor
from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver
from pyfea.solver.femm.domains.thermostatic.solver import FEMMThermostaticSolver

# Defines configuration file and solver output path
path_lib = "motors/prototype_1/optimization/configuration.uiv"
solver_folder = "motors/prototype_1/optimization/outputs"

# Defines the tubular linear motor via configuration parameters
TubularMotor = TubularLinearMotor(path_lib)

# Initializes FEA solvers (magnetostatic and thermostatic)
Magnetic = FEMMMagnetostaticSolver(solver_folder)
Thermal = FEMMThermostaticSolver(solver_folder)

static_results = static_evaluation(TubularMotor, Thermal, Magnetic)
print(static_results)