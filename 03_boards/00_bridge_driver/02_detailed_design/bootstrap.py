"""
Filename: bootstrap.py

Description:
    Calculates the size of the bootstrap capacitor 
    and the bootstrap resistor for the a single half bridge.
"""

from picounits import time, nullset, voltage, charge, current, nano, micro, kilo

# System Parameters
switch_frequency = 20 * kilo * (1/time)
duty_cycle = 0.5 * nullset

boot_current = 5 * current
resistor_drop = 0.2 * voltage
diode_forward = 1 * voltage

# Driver, Diode & Mosfet parameters
driver_voltage = 16 * voltage
driver_leakage = 50 * micro * current
driver_quiescent = 55 * micro * current
driver_minimum_voltage = 9 * voltage

max_gate_charge = 33 * nano * charge

# Computes the maximum voltage delta
delta_v = driver_voltage - diode_forward - driver_minimum_voltage - resistor_drop

# Computes the minimum capacitance at minimum frequency
frequency = 1 * (1 /time)
step_size = 500 * (1 / time)

# Calculates the minimum capacitance @ that switching frequency
Q_quiescent = driver_quiescent * duty_cycle / frequency
Q_leakage = driver_leakage * duty_cycle / frequency

# Total charge and required capacitor size
total = max_gate_charge + Q_leakage + Q_quiescent
boot_capacitor = total / delta_v

boot_resistor = resistor_drop / (boot_current)

print(f"Switch: {frequency:.3f}, Capacitor > {boot_capacitor:.3f}, Resistor: {boot_resistor:.3f}")

# Calculates the circuit time constant
time_constant = boot_resistor * boot_capacitor / duty_cycle
print(f"Boot circuit time constant: {time_constant:.3f}")
