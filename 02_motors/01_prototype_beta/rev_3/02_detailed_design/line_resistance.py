"""
Filename: line_resistance.py

Description:
    Calculates the wire cross-sectional
    area and ultimately the resistance
    for the motor power lines.
"""

from math import pi

from picounits import resolve_derived
from picounits import RESISTANCE, LENGTH, MILLI


# Loads the unit system
resolve_derived()

# Parameters
resistivity = 1.68 * 10 ** -8 * (RESISTANCE * LENGTH)

line_carrier = (10 * MILLI * LENGTH) * (10 * MILLI * LENGTH)
line_length = 500 * MILLI * LENGTH

wire_step = 0.1 * MILLI * LENGTH
wire_range = [0.5, 5] * MILLI * LENGTH

# Computed resistances & cross sectional area
steps = (wire_range[1] - wire_range[0]) / wire_step

for index in range(0, int(steps)+1):
    radius = wire_range[0] + index * wire_step

    # Computes the area and line carrier number
    area = pi * radius ** 2
    carrier_number = line_carrier / area

    # Computes the line resistance
    resistance = resistivity * line_length / area
    print(f"Radius: {radius}, Area: {area}, carrier_num: {carrier_number}, resistance: {resistance}")
