"""
Filename: velocity_impedance.py

Descriptions:
    Calculates the impedance of the motor per phase at
    a specific velocity and than calculates the maximal
    force that can be supplied.
    
    NOTE: Assumes low back-emf due to low frequency.
"""

from math import pi, sqrt

from picounits.core import unit_validator, quantities as q
from picounits.constants import (
    FORCE, CURRENT, IMPEDANCE, VOLTAGE, VELOCITY, MILLI, LENGTH, TIME, INDUCTANCE,
    CONDUCTIVITY
)

# Needs to be reworked for new model and also line to line is a better measure for max current.
# # Parameters
force_constant = 3.310 * FORCE / CURRENT
base_impedance = 4.574 * IMPEDANCE
pole_pitch = 20 * MILLI * LENGTH
max_velocity = 10 * VELOCITY

inductance = 32.014 * MILLI * INDUCTANCE
supply = 24 * VOLTAGE

slot_conductivity = 5.8 * 10 ** 7 * CONDUCTIVITY

@unit_validator(IMPEDANCE)
def line_resistance(wire_diameter: q, length: q) -> q:
    """ Calculates the line resistance of the wires """
    area = pi * (wire_diameter/2)**2 
    return length / (area * slot_conductivity)

@unit_validator(1/TIME)
def synchronous_freq(pitch: q, velocity: q) -> q:
    """ Calculates the synchronous frequency of the motor """
    return velocity / (2 * pitch)


@unit_validator(IMPEDANCE)
def inductive_reactance(frequency: q, inductance: q) -> q:
    """ Calculates the inductive reactance """
    return 2 * pi * frequency * inductance

base_impedance += line_resistance(1.02 * MILLI * LENGTH, 1 * LENGTH)

velocity_series = []
force_series = []
impedance_series = []
power_phase_series = []
power_total_series = []

for index in range(0, int(10 * max_velocity.value)):
    velocity = 0.1 * index * VELOCITY
    frequency = synchronous_freq(pole_pitch, velocity)
    inductive = inductive_reactance(frequency, inductance)

    impedance = (base_impedance**2 + inductive**2) ** 0.5
    current_phase = supply / impedance

    current_line = sqrt(3) * current_phase
    force = current_phase * force_constant
    power_phase = (current_phase**2) * base_impedance
    power_total = 3 * power_phase

    print(f"V: {velocity:.3f}, Z: {impedance:.3f}, I_phase: {current_phase:.3f}, "
          f"I_line: {current_line:.3f}, F: {force:.3f}, P_phase: {power_phase:.3f}, P_total: {power_total:.3f}")

    # Save series
    velocity_series.append(velocity.stripped)
    force_series.append(force.stripped)
    impedance_series.append(impedance.stripped)
    power_phase_series.append(power_phase.stripped)
    power_total_series.append(power_total.stripped)