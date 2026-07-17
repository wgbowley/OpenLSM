"""
Filename: main.py

Description:
    Reduced order model for tubular 
    linear synchronous motor. 
"""

from pathlib import Path
from picounits import Parser
from src.model import TubularMotor

# Loads unit system, material library & parameters
BASE_DIR = Path(__file__).parent
if not (BASE_DIR / "parameters.uiv").exists():
    raise FileNotFoundError("parameters.uiv not found in current directory")

materials = Parser.open(BASE_DIR / "lib/materials.uiv", BASE_DIR / "lib/metric.ut")
parameters = Parser.open(BASE_DIR / "parameters.uiv", BASE_DIR / "lib/metric.ut")

# Load/Constructs tubular linear motor model
model = TubularMotor(parameters, materials)
