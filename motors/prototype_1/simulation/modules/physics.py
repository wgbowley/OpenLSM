"""
Filename: physics.py

Description:
    Physics equation that describe motor
    dynamics such as clarke and park transforms,
    induced voltage, rk2_order_currents
"""

from math import pi, sqrt
from pyfea import Quantity as q, volt


def electrical_angle(displacement: q, pitch: q) -> q:
    """ 
    Calculates the electrical angle of the armature
    Ref. θ_e = π * displacement / pitch
    """
    return pi * displacement / pitch


def inverse_park_transform(currents: q, electrical_angle: q) -> q:
    """ converts d-p frame currents to stationary a-b frame """
    i_d, i_q = currents
    cos_t = electrical_angle.cos()
    sin_t = electrical_angle.sin()

    alpha = i_d * cos_t - i_q * sin_t
    beta  = i_d * sin_t + i_q * cos_t
    return (alpha, beta) * beta.unit


def inverse_clarke_transform(a_b_frame: q) -> q:
    """ Converts a-b stationary frame currents to three-phase (a,b,c) """
    alpha, beta = a_b_frame

    phase_a = alpha
    phase_b = 0.5 * (sqrt(3) * beta - alpha)
    phase_c = 0.5 * (-sqrt(3) * beta - alpha)
    return (phase_a, phase_b, phase_c) * phase_c.unit


def clarke_transform(a: q, b: q, c: q) -> q:
    """ Converts a, b, c values into alpha-beta components """
    alpha = (2 / 3) * (a - 0.5 * b - 0.5 * c)
    beta = (1 / sqrt(3)) * (b - c)

    return (alpha, beta) * beta.unit


def park_transform(a_b_frame: q, theta: q) -> q:
    """ Converts alpha-beta frame to d-q frame values """
    alpha, beta = a_b_frame
    cos_t = theta.cos()
    sin_t = theta.sin()
    
    d =  alpha * cos_t + beta * sin_t
    q = -alpha * sin_t + beta * cos_t
    return (d, q) * q.unit


def induced_voltage(delta_flux_linkage: q, delta_time: q) -> float:
    """
    Calculates total induced (back-EMF) voltage.
    E = Δψ / Δt
    """
    return delta_flux_linkage / delta_time


def differential_d_current(
    current_d: q, voltage_d: q, resistance: q, inductance: q, e_induced: q 
) -> q:
    """
    Differential equation for d-axis current.
    di_d/dt = (v_d - R*i_d - Δψ / Δt) / L (Assume uncoupled from q-axis)
    """
    return (voltage_d - resistance * current_d - e_induced) / inductance


def differential_q_current(
    current_q: q, voltage_q: q, resistance: q, inductance: q, e_induced: q 
) -> q:
    """
    Differential equation for q-axis current.
    di_q/dt = (v - R*i - Δψ / Δt) / L (Assume uncoupled from d-axis)
    """
    return (voltage_q - resistance * current_q - e_induced) / inductance


def rk_2nd_order_currents(
    currents: q, voltages: q, resistance: q, 
    inductance: q, induced: q, step_size: q
) -> q:
    """
    Solves the differential equations for the d-axis 
    and q-axis currents using ralston's method
    """
    # Extracts parameters from vectors
    i_d, i_q = currents[0], currents[1]
    v_d, v_q = voltages[0], voltages[1]
    h_d, h_q = inductance[0], inductance[1]
    # f_d, f_q = flux_linkages[0], flux_linkages[1]

    # Calculates the d-q axis induced voltages
    d_induced = 0 * volt
    q_induced = induced

    # Calculates the first step
    k1_d = differential_d_current(i_d, v_d, resistance, h_d, d_induced)
    k1_q = differential_q_current(i_q, v_q, resistance, h_q, q_induced)
    
    # Calculates the second step via time stepping
    i_ds = i_d + 3/ 4 * step_size * k1_d
    i_qs = i_q + 3/ 4 * step_size * k1_q
    
    k2_d = differential_d_current(i_ds, v_d, resistance, h_d, d_induced)
    k2_q = differential_q_current(i_qs, v_q, resistance, h_q, q_induced)
    
    # Final update using weighted average
    i_d += (1 / 3 * k1_d + 2 / 3 * k2_d) * step_size
    i_q += (1 / 3 * k1_q + 2 / 3 * k2_q) * step_size
    
    return (i_d, i_q) * i_q.unit