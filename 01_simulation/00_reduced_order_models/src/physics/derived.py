"""
Filename: derived.py

Description:
    Equations for derived parameters 
    within the reduced order model of 
    the motor.
"""

from math import pi, floor
from picounits import q


def compute_turns(length: q, thickness: q, wire_diameter: q, fill_factor: q) -> q:
    """ 
    Computes the number of turns while according for the insulation & stacking 
    Assumptions: Fill factor accounts for insulation & stacking
    """
    slot_section = length * thickness
    wire_section = pi * (wire_diameter / 2) ** 2

    effective_area = slot_section * fill_factor
    return floor(effective_area / wire_section)


def compute_inductance(turns: q, coil_len: q, mean_radius: q, permeability: q) -> q:
    """ Calculates the coils self-inductance independent of mutual inductance between coils """
    area = pi * mean_radius ** 2
    return (turns ** 2 * permeability * area) / coil_len


def compute_resistance(turns: q, mean_rad: q, wire_dia: q, resistivity: q) -> q:
    """ 
    Computes the slots resistance by using mean radius & conductor cross section 
    Assumptions: Uses mean radius as an approximate for approximate turn length
    """
    turn_length = 2 * pi * turns * mean_rad
    cross_section = pi * (wire_dia / 2) ** 2

    return resistivity * turn_length / cross_section


def compute_slot_volume(turns: q, wire_dia: q, mean_radius: q) -> q:
    """ 
    Computes the volume of the slots using wire diameter & mean radius 
    Assumptions: Uses mean radius as an approximate for approximate turn length
    """
    cross_section = pi * (wire_dia / 2) ** 2
    total_area = turns * cross_section
    mean_turn_length = 2 * pi * mean_radius

    return total_area * mean_turn_length
