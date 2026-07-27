"""
Filename: bootstrap.py

Description:
    Calculates the size of the bootstrap 
    capacitor and the bootstrap resistor 
    for the a single half bridge.
    
    NOTE: This model really needs to be updated
"""

from picounits import voltage, charge, micro, nano, time, kilo

# Parameters
f_sw = 20 * kilo * (1/time)
f_motor_start = 5 * (1/time)
duty_cycle = 0.5

driver_supply_voltage = 16 * voltage
diode_forward_voltage = 1 * voltage

# From the SUD90330E data-sheet
total_fet_gate_charge = 33 * nano * charge

# From the IRS2104 data-sheet
I_hbs = 50 * micro * (charge/time)
I_hb = 55 * micro * (charge/time)

HB_falling = 9 * voltage

# Computes the minimum capacitance for the system
delta_v_hb = driver_supply_voltage - diode_forward_voltage - HB_falling
Q_total = total_fet_gate_charge + (I_hbs * duty_cycle / f_motor_start) + (I_hb / f_motor_start)

boot_cap = Q_total / delta_v_hb

# Computes the boot resistor size
t_low_side_on = (1.0 - duty_cycle) / f_sw
boot_res = t_low_side_on / (5 * boot_cap)

energy_dissipated = 0.5 * boot_cap * (driver_supply_voltage - diode_forward_voltage) ** 2
current_peak = (driver_supply_voltage - diode_forward_voltage) / boot_res

print(f"Boot Capacitance: {boot_cap:.3f}, Boot Resistance: {boot_res:.3f}")
print(f"Energy Dissipated: {energy_dissipated:.3f}, I_peak: {current_peak:.3f}")
