"""
Filename: crystal_load_capacitors.py

Description:
    Calculates the capacitance for C1 and C2
    for the 8MHZ crystal (ECS-80-18-23B-JTN-TR)

    Reference:
    https://support.microchip.com/s/article/Calculating-crystal-load-capacitor
"""

from picounits.constants import CAPACITANCE, PICO

LOAD_CAP = 20 * PICO * CAPACITANCE
STRAY_CAP = 5 * PICO * CAPACITANCE   # Microchip.com assumes 5pf for this.


print(2 * (LOAD_CAP - STRAY_CAP))
