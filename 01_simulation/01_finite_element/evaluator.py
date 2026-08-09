"""
Filename: evaluator.py

Description:
    Evaluator for linear motors using PyFEA for magnetics
    parameter sweeps via finite element magnetic methods.
"""

from pathlib import Path
from pyfea.domain.units import Parser

from pcb_armature.model import PCBArmatureMotor

# Imports parameters from .uiv parameter file with units
BASE_DIR = Path(__file__).parents[0]
parameters_file = BASE_DIR / "pcb_armature/parameters.uiv"

# Imports the parameters (value:unit) into model
parameters = Parser.open(parameters_file)
motor = PCBArmatureMotor(parameters)
