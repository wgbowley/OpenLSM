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


def compute_resistance(turns: f, mean_rad: f, wire_dia: f, resistivity: f) -> f:
    """ 
    Computes the slots resistance by using mean radius & conductor cross section 
    Assumptions: Uses mean radius as an approximate for approximate turn length
    """
    turn_length = pi * turns * mean_rad
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