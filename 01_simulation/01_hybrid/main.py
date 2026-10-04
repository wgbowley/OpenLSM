"""
Filename: main.py

Description:
    Hybrid simulation for a tubular 
    linear synchronous motor.
"""

from pathlib import Path

from ifemm import Parser as iParser
from picounits import Parser

# Imports the parser and parses the .ans file
ROOT_DIR = Path(__file__).resolve().parents[0]
data = iParser.open(ROOT_DIR / 'resources/model.ans')

# Materials & Parameter files
parameters_path = ROOT_DIR / "parameters.uiv"
parameters = Parser.open(parameters_path, ROOT_DIR / "../derived.ut")

# (Work In Progress).
parameters.info()