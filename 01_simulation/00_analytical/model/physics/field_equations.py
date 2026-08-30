"""
Filename: field.py

Description:
    1D field equations for a analytical tubular linear motor.
"""

from builtins import float as f
from math import cosh


def _anti_derivative_function(z: f, pos: f, radius: f) -> f:
    """Computes the anti-derivative of the axial Biot-Savart pole field."""
    return (z - pos) / (radius ** 2 + (z - pos) ** 2) ** 0.5


def pole_field_strength(
    pos: f, translation: f, current: f, turns: f, length: f, radius: f
) -> f:
    """ Computes the field strength at a position z along the pole using the finite pole model """
    half_length = length / 2

    # Integral Bounds
    z_upper = translation + half_length
    z_lower = translation - half_length

    # Calculates the axial field components via integration
    term2 = _anti_derivative_function(z_upper, pos, radius)
    term1 = _anti_derivative_function(z_lower, pos, radius)

    h_term = turns * current / (2 * length)
    return h_term * (term2 - term1)


def _sech(z: float) -> float:
    """ Hyperbolic Secant Function """
    try:
        return 1 / cosh(z)

    except OverflowError:
        return 0.0


def _derived_pole_function(pos: f, offset: f, length: f, n: int) -> f:
    """ Computes the field strength at a position z along the pole """
    return _sech(n * (pos - offset) / length) ** 2


def dipole_field_strength(pos: f, start, h_field: f, length: f, n: int = 4) -> f:
    """ Computes the field strength at a position z along the dipole using a sech approximation. """
    half_length = length / 2

    # Pole locations
    z_upper = start - half_length
    z_lower = start + half_length

    # Calculates the axial field components (term 1 & term 2)
    term2 = _derived_pole_function(pos, z_upper, length, n)
    term1 = _derived_pole_function(pos, z_lower, length, n)

    # Calculates the field strength at `z_pos`
    return h_field * (term2 - term1)
