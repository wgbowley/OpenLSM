"""
Filename: field.py

Description:
    Equations for derived parameters 
    within the reduced order model of 
    the motor.
"""

from builtins import float as f
from math import pi, cos, sin, sqrt, floor


def electrical_angle(displacement: f, pitch: f) -> f:
    """ Calculates the electrical angle of the armature. """
    return pi * displacement / pitch


def inverse_park_transform(i_d: f, i_q: f, electrical_angle: f) -> tuple[f, f]:
    """ converts d-p frame currents to stationary a-b frame. """
    cos_t = cos(electrical_angle)
    sin_t = sin(electrical_angle)

    alpha = i_d * cos_t - i_q * sin_t
    beta  = i_d * sin_t + i_q * cos_t
    return alpha, beta


def inverse_clarke_transform(alpha_beta: tuple[f, f]) -> tuple[f, f, f]:
    """ Converts a-b stationary frame currents to three-phase (a,b,c). """
    alpha, beta = alpha_beta

    phase_a = alpha
    phase_b = 0.5 * (sqrt(3) * beta - alpha)
    phase_c = 0.5 * (-sqrt(3) * beta - alpha)
    return phase_a, phase_b, phase_c


def compute_turns(length: f, thickness: f, wire_diameter: f, fill_factor: f) -> f:
    """ Computes the number of turns while according for the insulation & stacking. """
    slot_section = length * thickness
    wire_section = pi * (wire_diameter / 2) ** 2

    effective_area = slot_section * fill_factor
    return floor(effective_area / wire_section)


def compute_inductance(turns: f, coil_len: f, mean_radius: f, permeability: f) -> f:
    """ Calculates the coils self-inductance independent of mutual inductance between coils. """
    area = pi * mean_radius ** 2
    return (turns ** 2 * permeability * area) / coil_len


def compute_resistance(turns: f, mean_rad: f, wire_dia: f, resistivity: f) -> f:
    """ Computes the slots resistance by using mean radius & conductor cross section. """
    turn_length = 2 * pi * turns * mean_rad
    cross_section = pi * (wire_dia / 2) ** 2

    return resistivity * turn_length / cross_section
