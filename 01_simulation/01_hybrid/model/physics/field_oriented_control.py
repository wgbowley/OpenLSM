"""
Filename: field_oriented_control.py

Description:
    FOC equations for the analytical tubular linear motor.
"""

from builtins import float as f
from math import pi, cos, sin, sqrt


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
