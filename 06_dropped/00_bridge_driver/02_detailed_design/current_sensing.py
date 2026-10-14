"""
Filename: current_sensing.py

Description:
    Calculates the max peak current for each 
    current sensing line based on their resistor size
"""

from picounits import Q, resolve_derived, expects
from picounits import resistance, voltage, current, milli

# Lazy imports derived units
resolve_derived()

# System Parameters
max_voltage_swing = 1 * voltage
gain = 20 * (voltage / voltage)     # -> (Nullset really)

# These are all R_2512_6332Metric if you need a different one find it.
line1 = 10 * milli * resistance
line2 = 5 * milli * resistance
line3 = 1 * milli * resistance


@expects(current)
def max_current(line: Q) -> Q:
    """ Calculates max line current based on a parameters """
    return max_voltage_swing / (line * gain)

# Computes the max currents
l1, l2, l3 = max_current(line1), max_current(line2), max_current(line3)

print(f"Max Currents: <l1: {l1}, l2: {l2}, l3: {l3}>")
