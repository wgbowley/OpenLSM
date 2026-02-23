"""
Filename: dynamic_physics.py

Description:
    Physics equation that describe motor
    dynamics such as clarke and park transforms,
    induced voltage, rk2_order_currents
"""

from math import pi, sqrt
from module.sim_definitions import MotorState

from pyfea import Quantity as q, weber
from pyfea.models.tubular_linear_motor.main import TubularLinearMotor

def electrical_angle(displacement: q, pitch: q) -> q:
    """  Calculates the electrical angle of the armature; θ_e = πd / p """
    return pi * displacement / pitch


def inverse_park_transform(currents: q, electrical_angle: q) -> q:
    """ converts d-p frame currents to stationary a-b frame """
    i_d, i_q = currents
    cos_t = electrical_angle.cos()
    sin_t = electrical_angle.sin()

    alpha = i_d * cos_t - i_q * sin_t
    beta  = i_d * sin_t + i_q * cos_t
    return (alpha, beta) * beta.unit


def park_transform(a_b_frame: q, theta: q) -> q:
    """ Converts alpha-beta frame to d-q frame values """
    alpha, beta = a_b_frame
    cos_t = theta.cos()
    sin_t = theta.sin()
    
    d =  alpha * cos_t + beta * sin_t
    q = -alpha * sin_t + beta * cos_t
    return (d, q) * q.unit


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


def induced_voltage(delta_flux_linkage: q, delta_time: q) -> q:
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
    di_d/dt = (v_d - R*i_d - E_d) / L_d
    """
    return (voltage_d - resistance * current_d - e_induced) / inductance


def differential_q_current(
    current_q: q, voltage_q: q, resistance: q, inductance: q, e_induced: q 
) -> q:
    """
    Differential equation for q-axis current.
    di_q/dt = (v_q - R*i_q - E_q) / L_q
    """
    return (voltage_q - resistance * current_q - e_induced) / inductance


def rk_2nd_order_currents(
    motor: TubularLinearMotor,
    state: MotorState, 
    prev_dq_flux: q,
    step_size: q
) -> q:
    """
    Solves the differential equations for the d-axis 
    and q-axis currents using Ralston's method
    """
    # Extract current d-q flux linkage from total flux
    a_b_frame = clarke_transform(
        state.total_flux[0], 
        state.total_flux[1], 
        state.total_flux[2]
    )
    elec_angle = electrical_angle(state.displacement, motor.pole_pitch)
    current_dq_flux = park_transform(a_b_frame, elec_angle)
    
    # Calculate induced voltages from flux change
    delta_flux_d = current_dq_flux[0] - prev_dq_flux[0]
    delta_flux_q = current_dq_flux[1] - prev_dq_flux[1]
    
    e_induced_d = induced_voltage(delta_flux_d, step_size)
    e_induced_q = induced_voltage(delta_flux_q, step_size)

    state.dq_induced = [e_induced_d.value, e_induced_q.value] * e_induced_d.unit
    
    # Extract parameters from state
    i_d, i_q = state.dq_currents[0], state.dq_currents[1]
    v_d, v_q = state.dq_voltages[0], state.dq_voltages[1]
    L_d, L_q = state.inductance[0], state.inductance[1]
    R = state.resistance
    
    # Calculates the first step (k1)
    k1_d = differential_d_current(i_d, v_d, R, L_d, e_induced_d)
    k1_q = differential_q_current(i_q, v_q, R, L_q, e_induced_q)
    
    # Calculates the second step via time stepping (k2 at 3/4 point)
    i_ds = i_d + 3/4 * step_size * k1_d
    i_qs = i_q + 3/4 * step_size * k1_q
    
    k2_d = differential_d_current(i_ds, v_d, R, L_d, e_induced_d)
    k2_q = differential_q_current(i_qs, v_q, R, L_q, e_induced_q)
    
    # Final update using weighted average (Ralston's method)
    i_d_new = i_d + (1/3 * k1_d + 2/3 * k2_d) * step_size
    i_q_new = i_q + (1/3 * k1_q + 2/3 * k2_q) * step_size
    
    return [i_d_new.value, i_q_new.value] * i_q_new.unit, current_dq_flux