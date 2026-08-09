"""
Filename: evaluator.py

Description:
    Evaluator for linear motors using PyFEA for magnetics
    parameter sweeps via finite element magnetic methods.
"""

from pathlib import Path
from pyfea.domain.units import Parser

from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver
from pcb_armature.model import PCBArmatureMotor

# Imports parameters from .uiv parameter file with units
BASE_DIR = Path(__file__).parents[0]
parameters_file = BASE_DIR / "pcb_armature/parameters.uiv"
solver_folder = BASE_DIR / "outputs"

# Imports the parameters (value:unit) into model
parameters = Parser.open(parameters_file)
motor = PCBArmatureMotor(parameters)

magnetic = FEMMMagnetostaticSolver(solver_folder)
domain = motor.construct_domain(magnetic)
magnetic.setup(domain)

### Generates the FEMM file for now which than can be edited within the FEMM GUI