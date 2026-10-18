"""
Filename: skin_depth.py

Description:
    Calculates the skin depth within a 
    conductor at a specific operating
    frequency.
"""

from math import pi

from picounits import resolve_derived
from picounits import RESISTANCE, INDUCTANCE, LENGTH, FREQUENCY


# Loads the unit system
resolve_derived()

# System Constants (Conductivity of 30 MS/m)
resistivity = 3.33 * 10 ** - 8 * (RESISTANCE * LENGTH)
free_space_permeability = 4 * pi * 10 ** -7 * (INDUCTANCE / LENGTH)

# System Parameters
conductor_permeability = 1

frequency_step = 10 * FREQUENCY
frequency_range = [10, 100] * FREQUENCY

steps = (frequency_range[1] - frequency_range[0]) / frequency_step


# Computes the skip depth over a range of frequencies
for index in range(0, int(steps)+1):
    frequency = frequency_range[0] + index  * frequency_step

    # Computes the skip depth at a specific frequency
    omega = 2 * pi * frequency

    permeability = conductor_permeability * free_space_permeability
    skip_depth = (2 * resistivity / (omega * permeability)) ** (1/2)

    print(f"Skip depth: {skip_depth:.3f} @ {frequency:.3f}")
