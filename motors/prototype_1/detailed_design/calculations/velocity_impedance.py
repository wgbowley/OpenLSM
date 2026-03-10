"""
Filename: velocity_impedance.py

Descriptions:
    Calculates the impedance of the motor per phase at
    a specific velocity and than calculates the maximal
    force that can be supplied.
    
    NOTE: Assumes low back-emf due to low frequency.
"""

from math import pi

from picounits.core import unit_validator, quantities as q
from picounits.constants import (
    FORCE, CURRENT, IMPEDANCE, VOLTAGE, VELOCITY, MILLI, LENGTH, TIME, INDUCTANCE
)

# Parameters
force_constant = 1.616 * FORCE / CURRENT
base_impedance = 1.831 * IMPEDANCE
pole_pitch = 20 * MILLI * LENGTH
max_velocity = 10 * VELOCITY

inductance = 15.83 * MILLI * INDUCTANCE
supply = 24 * VOLTAGE


@unit_validator(1/TIME)
def synchronous_freq(pitch: q, velocity: q) -> q:
    """ Calculates the synchronous frequency of the motor """
    return velocity / (2 * pitch)


@unit_validator(IMPEDANCE)
def inductive_reactance(frequency: q, inductance: q) -> q:
    """ Calculates the inductive reactance """
    return 2 * pi * frequency * inductance


velocity_series = []
force_series = []
impedance_series = []
for index in range(0, 10 * max_velocity.value):
    velocity =  0.1 * index * VELOCITY
    frequency = synchronous_freq(pole_pitch, velocity)
    inductive = inductive_reactance(frequency, inductance)
    
    impedance = (base_impedance**2 + inductive**2) ** 0.5
    current = supply / impedance
    force = current * force_constant

    
    print(f"V: {velocity:.3f}, R: {impedance:.3f}, I: {current:.3f}, F: {force:.3f}")
    velocity_series.append(velocity.stripped)
    force_series.append(force.stripped)
    impedance_series.append(impedance.stripped)
    