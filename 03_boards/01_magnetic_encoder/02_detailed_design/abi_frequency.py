"""
Filename: abi_frequency

Description:
    Calculates the a, b and index input frequencies
    depending on armature velocity.
    
    NOTE:
    Fore more reference: AS531_linear_magnetic_encoder.pdf 
"""

from picounits import MILLI, LENGTH, TIME

# Parameters
sample_num = 10
linear_speed_max = 650 * MILLI * LENGTH / TIME
linear_speed_step = 50 * MILLI * LENGTH / TIME

pole_length = 2 * MILLI * LENGTH

a_b_pulses_ratio = 256
states_per_pulse = 4

# Iterates until the max speed is reached
speed = linear_speed_step
while speed < linear_speed_max + linear_speed_step:
    # Calculates the index & AB frequency
    index_rate = speed/pole_length
    ab_edges = a_b_pulses_ratio * states_per_pulse

    ab_rate = ab_edges * index_rate
    print("-" * 12, f" {speed:.3f} ", "-" * 12)
    print(f"Fundamental: {index_rate:.3f}, {ab_rate:.3f}")

    # Calculates the sampling frequencies and sampling distance
    index_sample = index_rate*sample_num
    ab_sample = ab_rate * sample_num

    # sampling_distance = speed / ab_sample

    print(f"Sampling: {index_sample:.3f}, {ab_sample:.3f}")

    # Updates velocity for next iteration
    speed += linear_speed_step
