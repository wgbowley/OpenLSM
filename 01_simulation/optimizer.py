"""
Filename: optimizer.py

Description:
    converging random search using 1D analytical 
    simulation to refine the number of poles,
    slot size fraction and slot thickness.
"""

from pathlib import Path
from picounits import Parser

# Loads unit system, material library & parameters
ROOT_DIR = Path(__file__).resolve().parents[0]

# Materials & Parameter files
parameters_path = ROOT_DIR / "parameters.uiv"
parameters = Parser.open(parameters_path, ROOT_DIR / "derived.ut")
