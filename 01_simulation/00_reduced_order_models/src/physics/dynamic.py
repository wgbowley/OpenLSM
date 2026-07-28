"""
Filename: motor.py

Description:
    Dynamics equations for modelling the
    reduced order model of the motor.
"""

from builtins import float as f
from math import pi, cosh, cos, sin, sqrt


def compute_slot_z_field_strength(
    z_pos: f, z_start: f, current: f, turns: f, slot_len: f, slot_rad: f
) -> f:
    """
    Computes the field strength at a position z along the slot using the 
    finite continuous solenoid model in axial-symmetric modelling (Z, R).
    
    'z_slot' is in reference to the start of the slot. Hence center is shifted slot_len / 2
    """
    half_length = slot_len / 2
    center = z_start + half_length

    # Calculates the axial field components (term 1 & term 2)
    denom1 = ((z_pos - center + half_length) ** 2 + slot_rad ** 2) ** 0.5
    term1 = (z_pos - center + half_length) / denom1

    denom2 = ((z_pos - center - half_length) ** 2 + slot_rad ** 2) ** 0.5
    term2 = (z_pos - center - half_length) / denom2

    # Calculates maximal field_strength and returns position dependent strength
    h_term = turns * current / (2 * pi * slot_len)
    return h_term * (term1 - term2)


def compute_stator_z_field_strength(z_pos: f, z_start, pole_len: f, h_field: f, n: int = 4) -> f:
    """ 
    Computes the field strength at a position z along the stator using a sech approximation.
    finite continuous dipole model in axial-symmetric modelling (Z, R).
    
    'z_pole' is in reference to the start of the pole. Hence center is shifted pole_len / 2
    """
    def _sech(x):
        """ Hyperbolic Secant Function """
        try:
            return 1 / cosh(x)
        except OverflowError:
            return 0.0

    half_length = pole_len / 2
    center = z_start + half_length

    # Calculates the axial field components (term 1 & term 2)
    term1 = _sech(n * ((z_pos - center) + half_length) / pole_len) ** 2
    term2 = _sech(n * ((z_pos - center) - half_length) / pole_len) ** 2

    # Calculates the field strength at `z_pos`
    return h_field * (term1 - term2)


def electrical_angle(displacement: f, pitch: f) -> f:
    """ 
    Calculates the electrical angle of the armature
    Ref. θ_e = π * displacement / pitch
    """
    return pi * displacement / pitch


def inverse_park_transform(i_d: f, i_q: f, electrical_angle: f) -> tuple[f, f]:
    """ converts d-p frame currents to stationary a-b frame """
    cos_t = cos(electrical_angle)
    sin_t = sin(electrical_angle)

    alpha = i_d * cos_t - i_q * sin_t
    beta  = i_d * sin_t + i_q * cos_t
    return alpha, beta


def inverse_clarke_transform(alpha_beta: tuple[f, f]) -> tuple[f, f, f]:
    """ Converts a-b stationary frame currents to three-phase (a,b,c) """
    alpha, beta = alpha_beta

    phase_a = alpha
    phase_b = 0.5 * (sqrt(3) * beta - alpha)
    phase_c= 0.5 * (-sqrt(3) * beta - alpha)
    return phase_a, phase_b, phase_c
