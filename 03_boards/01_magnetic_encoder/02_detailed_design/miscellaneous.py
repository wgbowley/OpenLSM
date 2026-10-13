
"""
Filename: miscellaneous.py

Description:
    Computes/records methods for miscellaneous
    derived values.
"""

from picounits import TEMPERATURE, CURRENT, VOLTAGE, MILLI

# Absolute Ambient temperature range
max_temperature_range = abs(125 * TEMPERATURE - (-40 * TEMPERATURE))
tol_temperature_range = abs(70 * TEMPERATURE - (-30 * TEMPERATURE))
print(f"Max Temperature Range: {max_temperature_range:.3f}")
print(f"Tol (+/- 10um) Temperature Range: {tol_temperature_range:.3f}")


# Power usage
current = 16 * MILLI * CURRENT
voltage = 3.3 * VOLTAGE
power_usage = current * voltage
print(f"Typ. power usage: {power_usage:.3f}")

current = 21 * MILLI * CURRENT
voltage = 3.6 * VOLTAGE
power_usage = current * voltage
print(f"Max power usage: {power_usage:.3f}")
