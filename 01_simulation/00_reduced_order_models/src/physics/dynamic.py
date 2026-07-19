"""
Filename: motor.py

Description:
    Dynamics equations for modelling the
    reduced order model of the motor.
"""

from builtins import float as f
from math import pi, cosh


def compute_inductor_voltage(
    supply_voltage: f, current: f, resistance: f, induced_voltage: f
) -> f:
    """ Computes the potential difference across the inductor """
    voltage_drop = current * resistance
    return supply_voltage - voltage_drop + induced_voltage


def compute_current(
    current: f, voltage: f, inductance: f, resistance: f, time_step: f
) -> tuple[f, f]:
    """ 
    Computes the current using 2nd oder ralston's method however assumes the current changes.
    Assumptions: Assumes induced voltage term does not change (semi-implicit integration)
    """
    k1 = voltage / inductance

    # Updates the voltage for the next predicted frame
    voltage = voltage - resistance * 3 / 4 * k1 * time_step
    k2 = voltage / inductance

    di_dt = 1/3 * k1 + 2/3 * k2
    current += di_dt * time_step
    return current, di_dt


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
