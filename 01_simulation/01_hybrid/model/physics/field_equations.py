"""
Filename: field_equations.py

Description:
    Parametric winding layer equation and
    Axisymmetric biot-savart equation.
"""

import numpy as np

def standard_helix(t: np.ndarray, radius: float, pitch: float) -> np.ndarray:
    """ Generates the parametric helix array """
    x_wire = radius * np.cos(t)
    y_wire = radius * np.sin(t)
    z_wire = pitch / (2 * np.pi) * t

    # Creates an array of [(x1, y1, z1), ... , (xn, yn, zn)]
    return np.stack([x_wire, y_wire, z_wire], axis=1)


def derivative_helix(t: np.ndarray, radius: float, pitch: float, dt: float) -> np.ndarray:
    """" Generates the derivative array of the parametric helix"""
    dx = -radius * np.sin(t) * dt
    dy =  radius * np.cos(t) * dt
    dz =  (pitch / (2 * np.pi)) * dt * np.ones_like(t)

    # Creates an array of [(dx1, dy1, dz1), ... , (dxn, dyn, dzn)]
    return np.stack([dx, dy, dz], axis=1)
