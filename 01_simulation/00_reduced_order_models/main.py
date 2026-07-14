"""
Filename: main.py

Description:
    Reduced order model for tubular 
    linear synchronous motors. 
    
    ** TBD Exact Parameter Extraction **
"""

# This is a test to understand how to produce the b-field

from src.physics.motor import compute_slot_z_field_strength
from builtins import float as f
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from math import sin, pi


@dataclass
class MotorDynamics:
    " Holds motor dynamic variables "
    slot_pitch: float
    slot_length: float
    slot_radius: float
    turns: int
    rms_current: float
    a_current: float
    b_current: float
    c_current: float


def test_3_phase(t: f, state: MotorDynamics) -> MotorDynamics:
    """ Generates 3 phase currents in respect to time """
    state.a_current = state.rms_current * sin(t)
    state.b_current = state.rms_current * sin(t + 2*pi/3)
    state.c_current = state.rms_current * sin(t + 4*pi/3)
    
    return state

def armature_field(z_pos: float, state: MotorDynamics) -> float:
    """ Calculates the armature field at a z_position"""
    h_field = 0
    
    # Armature slots are [A+, B-, C+, A-, B+, C-]
    currents = [
        state.a_current, -state.b_current, state.c_current,
        -state.a_current, state.b_current, -state.c_current
    ]
    
    for i in range(6):
        z_slot = i * state.slot_pitch
        
        h_field += compute_slot_z_field_strength(
            z_pos, 
            z_slot, 
            currents[i], 
            state.turns, 
            state.slot_length, 
            state.slot_radius
        )
    
    # B=uH -> Hence U=u_o * u_m -> u_m = 1 hence -> B=u_o * H
    return h_field * (4 * pi * 10 ** -7)
    
    
state = MotorDynamics(
    slot_pitch=0.01, 
    slot_length=0.01, 
    slot_radius=0.005,
    turns=200,
    rms_current=2,
    a_current=0.0,
    b_current=0.0,
    c_current=0.0
)

t = pi / 2
state = test_3_phase(t, state)

print(f"At t={t}:")
print(f"  a_current = {state.a_current:.3f} A")
print(f"  b_current = {state.b_current:.3f} A")
print(f"  c_current = {state.c_current:.3f} A")

# Generate z positions
z_positions = np.linspace(-0.05, 0.150, 1000)
field_values = [armature_field(z, state) for z in z_positions]

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(z_positions * 1000, field_values, 'b-', linewidth=2)
ax.set_xlabel('Z Position (mm)')
ax.set_ylabel('Armature Field density (T)')
ax.set_title(f'Armature Field at t={t:.3f} s')

# Mark slot positions
for i in range(6):
    z_slot = i * state.slot_pitch - state.slot_length / 2
    ax.axvline(x=z_slot * 1000, color='red', linestyle='--', alpha=0.5)
    ax.text(z_slot * 1000, ax.get_ylim()[1]*0.9, f'S{i+1}', 
            ha='center', fontsize=9)

plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()