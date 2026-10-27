"""
Filename: potting_size.py

Description:
    Calculates the size of the potting
    around the copper assuming all copper 
    losses move through the potting interface.
"""

from math import pi, log
from picounits import THERMAL_CONDUCTIVITY, ENERGY, TEMPERATURE, POWER, LENGTH, MASS, MILLI

# Parameters
thermal_conductivity = 1.25 * THERMAL_CONDUCTIVITY
density = 1950 * MASS / LENGTH ** 3
approx_capacity = 1500 * ENERGY / (MASS * TEMPERATURE)

first_order_loss = 40 * POWER
outer_radius = 20 * MILLI * LENGTH
inner_radius = 11 * MILLI * LENGTH
axial_length = 60 * MILLI * LENGTH

for i in range(0, 9):
    # Calculates the temperature at the boundary of the potting
    term = 2 * pi * thermal_conductivity * axial_length
    DT = - first_order_loss / term * log(outer_radius / inner_radius)

    # Calculates the energy stored in that potting per k
    volume = pi * (outer_radius ** 2 - inner_radius ** 2) * axial_length
    mass = volume * density

    capacity = approx_capacity * mass

    print(f"{outer_radius}, {inner_radius}, {volume}: {DT}, {capacity}")
    outer_radius -= 1 * MILLI * LENGTH
