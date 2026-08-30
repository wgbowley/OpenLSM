"""
Filename: ina214a_ref.py

Description:
    Calculates the v_out based on
    the components referencing.
"""

from pathlib import Path

from picounits import Parser
from picounits import Q, expects, NULLSET, CURRENT, VOLTAGE, IMPEDANCE, MILLI

# Imports derived units
ROOT_DIR = Path(__file__).resolve().parents[0]
parameters = Parser.import_derived(ROOT_DIR / "../../metric.ut")

# System parameters
current_size = 300 * MILLI * CURRENT
current_range = [0, 3] * CURRENT
shunt = 10 * MILLI * IMPEDANCE

ref = 2 * VOLTAGE
gain = 50 * NULLSET


@expects(VOLTAGE)
def output(gain: Q, v_d: Q, v_ref: Q) -> Q:
    """ Computes the output voltage of the ina241 based on parameters """
    return gain * v_d + v_ref / 2


# Calculates the loop numbers
steps = (current_range[1] - current_range[0]) // current_size

for index in range(0, int(steps)+1):
    # Computes the sampling current & shunt drop
    current = current_size * index
    voltage_drop = current * shunt

    # Computes and prints the respective output voltage
    output_voltage = output(gain, voltage_drop, ref)
    print(f"current: {current:.3f}, drop: {voltage_drop:.3f}, output: {output_voltage:.3f}")
