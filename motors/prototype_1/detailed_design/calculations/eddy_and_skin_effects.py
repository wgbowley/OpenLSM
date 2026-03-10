"""
Filename: eddy_and_skin_effects.py

Descriptions:
    Calculates the skin depth of materials
    at different synchronous frequency and
    the eddy currents within the material from 
    radial b-field data. 
    
    NOTE:
    The radial b-field data is field at phase angle pi/2 at 1A.
    Hence pa = 1, pb=-0.5, bc=-0.5
"""

from math import pi

from picounits.core import unit_validator, quantities as q
from picounits.constants import (
    VELOCITY, MILLI, LENGTH, TIME, MICRO, MASS, POWER,
    CONDUCTIVITY, PERMEABILITY, FLUX_DENSITY
)

FREE_SPACE = 4 * pi * 10 ** -7 * PERMEABILITY
EDDY_CO = 10 * MICRO * CONDUCTIVITY

# Parameters
pole_pitch = 20 * MILLI * LENGTH
max_velocity = 10 * VELOCITY

# Core is cf-nylon
core_permeability = 1 * FREE_SPACE
core_conductivity = 0.1 * CONDUCTIVITY
core_volume = 9.387 * MICRO * LENGTH ** 3
core_thickness = 4 * MILLI * LENGTH
core_eddy_co_efficient = 0 * EDDY_CO

# Slot is copper
slot_b_max = 0.8 * FLUX_DENSITY
slot_permeability = 1 * FREE_SPACE
slot_conductivity = 5.8 * 10 ** 7 * CONDUCTIVITY
total_slot_volume = 9.839 * MICRO * LENGTH ** 3
volume_per_slot = total_slot_volume / 12
slot_thickness = 3 * MILLI * LENGTH
slot_eddy_co_efficient = 6.67 * 10 ** 7 * EDDY_CO

# Heat Sink is aluminum
sink_b_max = 0.2 * FLUX_DENSITY
sink_permeability = 1 * FREE_SPACE
sink_conductivity = 3.5 * 10 ** 7 * CONDUCTIVITY
heat_sink_volume = 21.488 * MICRO * LENGTH ** 3
sink_thickness = 5 * MILLI * LENGTH
sink_eddy_co_efficient = 1.30 * 10 ** 7 * EDDY_CO


@unit_validator(1/TIME)
def synchronous_freq(pitch: q, velocity: q) -> q:
    """ Calculates the synchronous frequency of the motor """
    return velocity / (2 * pitch)


@unit_validator(LENGTH)
def skin_depth(conductivity: q, permeability: q, frequency: q) -> q:
    """ Calculates the skin depth of the component """
    res = 1/conductivity
    return (2 * res / (2*pi*frequency * permeability)) ** 0.5


@unit_validator(POWER)
def eddy_losses(ke: q, b_max: q, freq: q, thickness:q, volume: q) -> q:
    """ Calculates the eddy losses within different components """
    return ke * b_max ** 2 * freq ** 2 * thickness ** 2 * volume


for index in range(1, 10 * max_velocity.value):
    velocity = 0.1 * index * VELOCITY
    frequency = synchronous_freq(pole_pitch, velocity)
    
    # Calculates component skin depth and eddy losses at sync frequency
    sink_d = skin_depth(sink_conductivity, sink_permeability, frequency)
    sink_p = eddy_losses(
        sink_eddy_co_efficient, sink_b_max, frequency, sink_thickness, heat_sink_volume 
    )
    
    slot_d = skin_depth(slot_conductivity, slot_permeability, frequency)
    slot_p = eddy_losses(
        slot_eddy_co_efficient, slot_b_max, frequency, slot_thickness, volume_per_slot
    )
    
    print(f"Slot(V: {velocity:.3f}, F: {frequency:.3f}, D: {slot_d:.3f}, P: {slot_p:.3f}")
    print(f"HeatSink(V: {velocity:.3f}, F: {frequency:.3f}, D: {sink_d:.3f}, P: {sink_p:.3f}")
