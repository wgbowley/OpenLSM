"""
Filename: derived.py

Description:
    Equations for derived parameters 
    within the reduced order model of 
    the motor.
"""

from math import pi, floor
from builtins import float as f


def compute_turns(length: f, thickness: f, wire_diameter: f, fill_factor: f) -> f:
    """ 
    Computes the number of turns while according for the insulation & stacking 
    Assumptions: Fill factor accounts for insulation & stacking
    """
    slot_section = length * thickness
    wire_section = pi * (wire_diameter / 2) ** 2

    effective_area = slot_section * fill_factor
    return floor(effective_area / wire_section)


def compute_inductance(turns: f, coil_len: f, mean_radius: f, permeability: f) -> f:
    """ Calculates the coils self-inductance independent of mutual inductance between coils """
    area = pi * mean_radius ** 2
    return (turns ** 2 * permeability * area) / coil_len


def compute_resistance(turns: f, mean_rad: f, wire_dia: f, resistivity: f) -> f:
    """ 
    Computes the slots resistance by using mean radius & conductor cross section 
    Assumptions: Uses mean radius as an approximate for approximate turn length
    """
    turn_length = 2 * pi * turns * mean_rad
    cross_section = pi * (wire_dia / 2) ** 2

    return resistivity * turn_length / cross_section


def compute_slot_volume(turns: f, wire_dia: f, mean_radius: f) -> f:
    """ 
    Computes the volume of the slots using wire diameter & mean radius 
    Assumptions: Uses mean radius as an approximate for approximate turn length
    """
    cross_section = pi * (wire_dia / 2) ** 2
    total_area = turns * cross_section
    mean_turn_length = 2 * pi * mean_radius

    return total_area * mean_turn_length


def compute_stator_field_strength(flux_density: f, relative_permeability: f, permeability: f) -> f:
    """
    Computes the stator field strength using B=u_o * u_r * H relationship.
    Assumptions: 
    This model works well for ironless machines as flux density and field strength are proportional
    due to a lack of non-linear materials in them.
    """
    return flux_density / (relative_permeability * permeability)
