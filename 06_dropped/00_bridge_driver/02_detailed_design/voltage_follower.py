"""
Filename: voltage_follower.py

Description:
    Calculates a resistor divider for generating 
    a reference voltage from a given supply voltage.
"""

from pathlib import Path
from picounits import Parser, Q, expects

from picounits import VOLTAGE, RESISTANCE

# Imports derived units
ROOT_DIR = Path(__file__).resolve().parents[0]
parameters = Parser.import_derived(ROOT_DIR / "../../derived.ut")

# System Parameters
input_voltage = 3.3 * VOLTAGE
output_voltage = 2.0 * VOLTAGE

# Parameter Space
resistors = [
    10, 22, 47, 100, 150, 220, 330, 470,
    560, 680, 820, 1000, 1100, 1200, 1300, 1500,
    1600, 1800, 2000, 2200, 2400, 2700, 3000, 3300,
    3600, 3900, 4300, 4700, 5100, 5600, 6200, 6800,
    7500, 8200, 9100, 10000
] * RESISTANCE


@expects(RESISTANCE)
def divider_resistance(vin: Q, vout: Q, r_2: Q) -> Q:
    """ Calculates the resistance of r_1 """
    return (r_2 * vin - vout * r_2 ) / vout

for item in resistors:
    r1 = divider_resistance(input_voltage, output_voltage, item)
    print(f"r1: {r1}, r2: {item}")
