"""
Filename: potting_size.py

Description:
    Calculates the size of the heat sink
    around the potting move through the 
    potting-heat sink interface
"""

from math import pi, log
from picounits import THERMAL_CONDUCTIVITY, TEMPERATURE, POWER, LENGTH, MILLI

# Parameters
first_order_loss = 40 * POWER
inner_radius = 12.30 * MILLI * LENGTH
outer_radius = 25 * MILLI * LENGTH
axial_length = 60 * MILLI * LENGTH

k_aluminum = 167 * THERMAL_CONDUCTIVITY

# Temperature from potting_size.py script
temp_potting_boundary = -7.386 * TEMPERATURE


for i in range(0, 9):
    # Thermal resistance of the heat sink
    r_aluminum = log(outer_radius / inner_radius) / (2 * pi * k_aluminum * axial_length)
    dT_aluminum = - first_order_loss * r_aluminum

    # Calculates the drop temperature over copper to potting to 6061 interface.
    temperature_drop = temp_potting_boundary + dT_aluminum

    print(f"{outer_radius}, {inner_radius}: {dT_aluminum}, {temperature_drop}")
    outer_radius -= 1 * MILLI * LENGTH
