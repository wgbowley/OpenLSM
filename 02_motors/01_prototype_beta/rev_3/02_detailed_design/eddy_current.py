"""
Filename: eddy_current.py

Description:
    Calculates the standard depth of penetration
    and proportionality between variables.
"""

from math import pi

from picounits import Q, expects
from picounits import LENGTH, TIME, CONDUCTIVITY, PERMEABILITY, MEGA, MILLI


# (Note). This doesn't seem correct for the application

@expects(LENGTH)
def skin_depth(frequency: Q, u_r: Q, conductivity: Q) -> Q:
    """ Calculates the skin depth of penetration within the armature. """
    u_0 = 4 * pi * 10 ** - 7 * PERMEABILITY

    return 1 / (pi * frequency * u_r * u_0 * conductivity) ** 0.5


# Variables
velocity = 650 * MILLI * (LENGTH / TIME)
pole_pitch = 20 * MILLI * LENGTH

permeability = 1
armature_conductivity = 25 * MEGA * CONDUCTIVITY
electrical_frequency = velocity / pole_pitch

depth = skin_depth(electrical_frequency, permeability, armature_conductivity)

print(f"depth: {depth}, {2 * depth}, {3 * depth}")
