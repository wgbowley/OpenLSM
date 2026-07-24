"""
Filename: thermistor_basis.py

Description:
    Calculates the basis resistor for the 10kΩ thermistor 
    divider circuit into the Multiplexer (74HC4051) and 
    than ADC (PB0 STM32G431K8Tx)
"""

from math import floor, log

from picounits import Q
from picounits import IMPEDANCE, VOLTAGE, TIME, TEMPERATURE, CAPACITANCE, NANO, PICO

# Parameters
supply_voltage = 3.3 * VOLTAGE
target_frequency = 100 * (1/TIME)
adc_range = 4096

# Thermistor basis resistor
thermistor_basis = 1800 * IMPEDANCE

# Decoupling capacitor
decoupling = 10 * NANO * CAPACITANCE

# Multiplexer
mux_resistance = 100 * IMPEDANCE
mux_capacitance = 5 * PICO * CAPACITANCE

# 0c to 150c
OPERATING_RANGE = [32624, 181.370] * IMPEDANCE


def thermistor_range(basis_res: Q) -> tuple[Q, Q]:
    """ Returns the voltage range & ADC step range """
    values = []
    for res in OPERATING_RANGE:
        # Calculates the voltage output at that temperature
        res_network = basis_res / (res + basis_res)
        values.append(supply_voltage * res_network)

    # Calculates the voltage range and ADC step range
    voltage_range = values[-1]-values[0]
    step_range = voltage_range / (supply_voltage) * adc_range

    return voltage_range, floor(step_range)


def dc_impedance(basis_res: Q) -> tuple[Q, Q]:
    """ Computes the dc input impedance at (PB0 STM32G431K8Tx) """
    values = []
    for res in OPERATING_RANGE:
        # Calculates the thvenin resistance of the divider
        r_thvenin = (basis_res * res) / (basis_res + res)
        values.append(r_thvenin)

    return values


# Calculates results and displays results
voltage_range, step_range = thermistor_range(thermistor_basis)
temperature_impedance = dc_impedance(thermistor_basis)

print(f"Operating range: {273 * TEMPERATURE} to {(273 + 150) * TEMPERATURE}")
print(f"ADC voltage range: {voltage_range:.3f}")
print(f"ADC step range: {step_range:.3f}, Temp/Step: {150 * TEMPERATURE / step_range:.3f}")

# Calculates time constant
for res in (temperature_impedance):
    total_cap = decoupling + mux_capacitance
    tau = (res + mux_resistance) * total_cap

    print(f"Time required at {res:.3f}, {tau * log(4096):.3f}")

# Calculates thermal load
for res in OPERATING_RANGE:
    # Calculates the current through the divider
    current = supply_voltage / (res + thermistor_basis)

    print(f"Thermistor basis loss: {current ** 2 * thermistor_basis:.3f}")
    print(f"Thermistor Power loss: {current**2 * res:.3f}")
    print(f"Total Resistance power loss: {current**2 * (res + thermistor_basis):.3f}")
