"""
Filename: field_equations.py

Description:
    Parametric winding layer equation and
    Axisymmetric biot-savart equation.
"""

from builtins import float as f

import numpy as np


def _anti_derivative_function(z: f, pos: f, radius: f) -> f:
    """ Computes the anti-derivative of the axial Biot-Savart pole field. """
    return (z - pos) / (radius ** 2 + (z - pos) ** 2) ** 0.5


def _integrand_pole_model(pos: f, length: f, radius: f) -> f:
    """ Computes the integrand along the z-axis using the finite pole model """
    half_length = length / 2

    # Calculates the axial field components via integration
    term2 = _anti_derivative_function(half_length, pos, radius)
    term1 = _anti_derivative_function(-half_length, pos, radius)

    return term2 - term1


def kernel_limit(dropoff: f, length: f, radius: f) -> f:
    """ Calculates the limit for the kernel with dropoff being the limit """
    original = _integrand_pole_model(0.0, length, radius)
    target = dropoff * original

    # Restricts the search space within bounds [radius, inf]
    high = radius
    while _integrand_pole_model(high, length, radius) > target:
        high *= 2

    # Uses the bisection search to refine the result
    low = radius
    for _ in range(32):
        center = (low + high) / 2
        value = _integrand_pole_model(center, length, radius)

        if value > target:
            low = center
        else:
            high = center

    return (low + high) / 2


def standard_helix(t: np.ndarray, radius: float, length: float, turns: float) -> np.ndarray:
    """ Generates the parametric helix array """
    x_wire = radius * np.cos(t)
    y_wire = radius * np.sin(t)
    z_wire = length / (2 * np.pi * turns) * t

    # Creates an array of [(x1, y1, z1), ... , (xn, yn, zn)]
    return np.stack([x_wire, y_wire, z_wire], axis=1)


def derivative_helix(t: np.ndarray, radius: float, length: float, turns: float, dt: float) -> np.ndarray:
    """" Generates the derivative array of the parametric helix"""
    dx = -radius * np.sin(t) * dt
    dy =  radius * np.cos(t) * dt
    dz =  (length / (2 * np.pi * turns)) * dt * np.ones_like(t)

    # Creates an array of [(dx1, dy1, dz1), ... , (dxn, dyn, dzn)]
    return np.stack([dx, dy, dz], axis=1)


def biot_sum_integrand(r_eval: np.ndarray, r_wire: np.ndarray, dl: np.ndarray) -> np.ndarray:
    """ 3D Biot savart sum of the integrand """
    r_vec = r_eval[..., None, :] - r_wire

    # Calculates magnitude of r_vec and forces safe flooring
    distance = np.linalg.norm(r_vec, axis=-1)
    flooring = np.maximum(distance, 1e-10)

    # Cross and resulting summation of the integrand
    cross = np.linalg.cross(dl, r_vec)
    return np.sum(cross / flooring[..., None]**3, axis=-2)
