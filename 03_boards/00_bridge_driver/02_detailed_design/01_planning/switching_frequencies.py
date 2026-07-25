"""
Filename: switching_frequencies.py

Description:
    Calculates the maximum switching frequency based on the 
    fet using total gate charge and driver max supply current.
    
    NOTE:
    Assumes no parasitic elements are affecting the raise/fall times. 
"""

from picounits.constants import charge, time, milli, nano


# parameters (max. not typ. or min.)
driver_turn_on_delay = 820 * nano * time
driver_turn_off_delay = 220 * nano * time

driver_turn_on_raise = 170 * nano * time
driver_turn_off_fall = 90 * nano * time

# typ. high side pulse current
driver_peak_current = 210 * milli * (charge / time)

fet_turn_on_delay = 5.2 * nano * time
fet_turn_off_delay = 17 * nano * time
fet_turn_on_raise = 3 * nano * time
fet_turn_off_fall = 2.5 * nano * time

fet_total_gate_charge = 25 * nano * charge

# Computes the total switching time and than the maximum frequency of the system
driver = driver_turn_on_delay + driver_turn_off_delay + driver_turn_on_raise + driver_turn_off_fall
fet = fet_turn_on_delay + fet_turn_off_delay + fet_turn_on_raise + fet_turn_off_fall
network_time = driver + fet

fet_charge_time = fet_total_gate_charge / driver_peak_current 
print(f"Maximum possible charge time: {fet_charge_time:.3f}")

total_time = network_time + fet_charge_time
print(f"Maximum possible frequency: {1/total_time:.3f} @ {driver_peak_current:.3f}")
