"""
Filename: dynamic_physics.py

Description:
    Physics equation that describe motor
    dynamics such as clarke and park transforms,
    induced voltage, rk2_order_currents
"""

from math import pi, sqrt
from module.sim_definitions import MotorState

from pyfea import Quantity as q
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


def _id_dt(ms: MotorState, d_v: q, d_i: q, q_i: q, f_e: q) -> q:
    """ Differential equation for d-axis current with q coupling """
    standard = d_v - ms.resistance * d_i
    coupling = f_e * q_i * ms.inductance[1] 

    return (standard + coupling) / ms.inductance[0]


def _id_qt(ms: MotorState, q_v: q, q_i: q, d_i: q, f_e: q, f_m: q) -> q:
    """ Differential equation for q-axis current with d coupling """
    standard = q_v - ms.resistance * q_i
    coupling = f_e * (d_i * ms.inductance[0] + f_m)
    
    return (standard - coupling) / ms.inductance[1]


def rk_2nd_order_currents(
    motor: TubularLinearMotor, state: MotorState, step_size: q, magnet_flux: q
) -> None:
    """ 
    Solves the differential equations for the d-axis and q-axis 
    currents using ralston's method
    """
    ms, f_m = state, magnet_flux
    f_e = pi * ms.velocity / motor.pole_pitch
    
    # Calculates the first & second step
    k1_d = _id_dt(ms, ms.dq_voltages[0], ms.dq_currents[0], ms.dq_currents[1], f_e)
    k1_q = _id_qt(ms, ms.dq_voltages[1], ms.dq_currents[1], ms.dq_currents[0], f_e, f_m)

    # Calculates the second step
    new_velocity = ms.velocity + 3 / 4 * ms.acceleration * step_size
    f_e = pi * new_velocity / motor.pole_pitch
 
    i_ds = ms.dq_currents[0] + 3 / 4 * step_size * k1_d
    i_qs = ms.dq_currents[1] + 3 / 4 * step_size * k1_q

    k2_d = _id_dt(ms, ms.dq_voltages[0], i_ds, i_qs, f_e)
    k2_q = _id_qt(ms, ms.dq_voltages[1], i_qs, i_ds, f_e, f_m)

    # Updates current with weighted average derivative
    ms.dq_currents[0] += (1 / 3 * k1_d + 2 / 3 * k2_d) * step_size
    ms.dq_currents[1] += (1 / 3 * k1_q + 2 / 3 * k2_q) * step_size
 