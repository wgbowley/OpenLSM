"""
Filename: fundamental_frequencies.py

Description:
    This script generates a table of linear velocity to fundamental 
    frequencies for the Prototype Alpha linear motor.
"""

from picounits import VELOCITY, LENGTH, MILLI


# Mechanical Parameters
pole_to_pole = 14 * MILLI * LENGTH
pole_pitch = 2 * pole_to_pole

# Control Parameters
velocity_step = 100 * MILLI * VELOCITY
velocity_max = 1000 * MILLI * VELOCITY

velocity = velocity_step
while velocity < velocity_max:
    frequency = velocity / pole_pitch
    
    print(f"Velocity: {velocity:.3f}, frequency: {frequency:.3f}")
    velocity += velocity_step
