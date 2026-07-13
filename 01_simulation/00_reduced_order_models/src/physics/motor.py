"""
Filename: motor.py

Description:
    Dynamics equations for modelling the
    reduced order model of the motor.
"""

from builtins import float as f
from math import pi


def compute_inductance(turns: f, coil_len: f, mean_radius: f, permeability: f) -> f:
    """ Calculates the coils self-inductance independent of mutual inductance between coils """
    area = pi * mean_radius ** 2
    return (turns ** 2 * permeability * area) / coil_len


def computes_inductor_voltage(
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
    
    di_dt = (1/3 * k1 + 2/3 * k2)
    current += di_dt * time_step
    return current, di_dt   
    

def compute_z_field_strength(
    z_pos: f, z_slot: f, current: f, turns: f, slot_len: f, slot_rad: f
) -> f:
    """
    Computes the field strength at a position z within or outside of the slot using the 
    finite continuous solenoid model in axial-symmetric modelling (Z, R).
    
    'z_slot' is in reference to the front of the slot. Hence center is shifted back slot_len / 2
    """
    half_length = slot_len / 2
    slot_center = z_pos - (z_slot - half_length)
    
    # Calculates the axial field components (term 1 & term 2)
    denom1 = slot_len * (slot_rad ** 2 + (slot_center + half_length)**2) ** 0.5
    term1 = (slot_center + half_length) / denom1
    
    denom2 = slot_len * (slot_rad ** 2 + (slot_center - half_length)**2) ** 0.5
    term2 = (slot_center - half_length) / denom2
    
    # Calculates maximal field_strength and returns position dependent strength
    h_term = turns * current / 2 
    return h_term * (term1 - term2)