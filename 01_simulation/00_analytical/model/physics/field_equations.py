"""
Filename: field.py

Description:
    Field Equations for modelling 
    the tubular linear motor.
"""

from builtins import float as f
from math import cosh, pi


def compute_pole_field_strength(
    pos: f, translation: f, current: f, turns: f, length: f, radius: f
) -> f:
    """
    Computes the field strength at a position z along the pole using
    the finite pole model in axial-symmetric modelling (Z-R)
    """
    half_length = length / 2

    # Integral Bounds
    z_upper = translation + half_length
    z_lower = translation - half_length

    # Calculates the axial field components via integration
    term1 = (pos + z_upper) / (radius ** 2 + (pos + z_upper) ** 2) ** 0.5
    term2 = (pos + z_lower) / (radius ** 2 + (pos + z_lower) ** 2) ** 0.5

    # Calculates maximal field strength and returns position dependent strength
    # Need to derive the formula again to check if pi cancels out or not.
    h_term = turns * current / (2 * pi * length)
    return h_term * (term2 - term1)


def compute_dipole_field_strength(
    pos: f, start, h_field: f, length: f, n: int = 4
) -> f:
    """ 
    Computes the field strength at a position z along the dipole using a sech approximation.
    """
    def _sech(x):
        """ Hyperbolic Secant Function """
        try:
            return 1 / cosh(x)
        except OverflowError:
            return 0.0

    half_length = length / 2

    # Calculates the axial field components (term 1 & term 2)
    term1 = _sech(n * ((pos - start) + half_length) / length) ** 2
    term2 = _sech(n * ((pos - start) - half_length) / length) ** 2

    # Calculates the field strength at `z_pos`
    return h_field * (term1 - term2)
