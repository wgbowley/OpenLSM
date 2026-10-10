"""
Filename: main.py

Description:
    Hybrid simulation for a tubular 
    linear synchronous motor.
"""

from math import sin, pi
from pathlib import Path

from ifemm import Parser as iParser
from picounits import Parser

from model.solver import Solver

# Imports the parser and parses the .ans file
ROOT_DIR = Path(__file__).resolve().parents[0]
data = iParser.open(ROOT_DIR / 'resources/model.ans')

# Materials & Parameter files
parameters_path = ROOT_DIR / "parameters.uiv"
parameters = Parser.open(parameters_path, ROOT_DIR / "../derived.ut")

solver = Solver(parameters, data)

solver.i_pha, solver.i_phb, solver.i_phc = 10 * sin(0), 10 * sin(2 * pi / 3), 10 * sin(4 * pi / 3)

