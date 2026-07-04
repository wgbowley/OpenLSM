"""
Filename: thermistor_basis_resistor.py

Description:
    Calculates the basis resistor for the 10kOhm thermistor
    divider circuit into the 74HC4051 and than ADC
"""

from picounits.core import quantities as q, unit_validator
from picounits.constants import (
    IMPEDANCE, VOLTAGE, CURRENT, MILLI, NANO, CAPACITANCE, TIME, PICO
)

SUPPLY_VOLTAGE = 3.3 * VOLTAGE
ADC_RANGE = 4096

MAXIMUM_SELF_HEAT = 100 * MILLI * (VOLTAGE * CURRENT)

DECODING_CAP = 10 * NANO * CAPACITANCE
FREQUENCY = 100 * (1 / TIME)

# IDEAL_IMPEDANCE = 1000 * IMPEDANCE; Generic rule. Not physics
# 3 time constants to max; Half the allowed time to reach. hence "6 steps"
IDEAL_TIME_CONSTANT = 1/(6*(FREQUENCY))

MUX_RESISTANCE = 100 * IMPEDANCE
MUX_CAPACITANCE = 5 * PICO * CAPACITANCE


THERMISTOR_LOOKUP_TABLE = [
    664169, # -50
    334274, # -40
    176133, # -30
    96761,  # -20
    55218,  # -10
    32624,  # 0
    19897,  # 10
    12493,  # 20
    8056.0, # 30
    5323.9, # 40
    3598.7, # 50
    2483.8, # 60
    1747.7, # 70
    1251.8, # 80
    911.59, # 90
    674.11, # 100
    505.68, # 110
    384.41, # 120
    295.881, # 130
    230.400, # 140
    181.370, # 150
] * IMPEDANCE

OPERATING_RANGE = [19897, 181.370] * IMPEDANCE # 10 -> 150c

standard_resistor_series = [
    1, 10, 100, 1.1, 11, 110, 1.2, 12, 120, 1.3, 13, 130, 1.5, 15, 150, 
    1.6, 16, 160, 1.8, 18, 180, 2, 20, 200, 2.2, 22, 220, 2.4, 24, 240, 
    2.7, 27, 270, 3, 30, 300, 3.3, 33, 330, 3.6, 36, 360, 3.9, 39, 390, 
    4.3, 43, 430, 4.7, 47, 470, 5.1, 51, 510, 5.6, 56, 560, 6.2, 62, 620, 
    6.8, 68, 680, 7.5, 75, 750, 8.2, 82, 820, 9.1, 91, 910, 1000, 10000, 
    100000, 1100, 11000, 110000, 1200, 12000, 120000, 1300, 13000, 130000, 
    1500, 15000, 150000, 1600, 16000, 160000, 1800, 18000, 180000, 2000, 
    20000, 200000, 2200, 22000, 220000, 2400, 24000, 240000, 2700, 27000, 
    270000, 3000, 30000, 300000, 3300, 33000, 330000, 3600, 36000, 360000, 
    3900, 39000, 390000, 4300, 43000, 430000, 4700, 47000, 470000, 5100, 
    51000, 510000, 5600, 56000, 560000, 6200, 62000, 620000, 6800, 68000, 
    680000, 7500, 75000, 750000, 8200, 82000, 820000, 9100, 91000, 910000, 
    1000000, 3000000, 9100000, 1100000, 3300000, 10000000, 1200000, 3600000, 
    12000000, 1300000, 3900000, 13000000, 1500000, 4300000, 15000000, 1600000, 
    5100000, 16000000, 1800000, 5600000, 18000000, 2000000, 6200000, 20000000, 
    2200000, 6800000, 22000000, 2400000, 7500000, 2700000, 8200000, 0.5, 0.22
] * IMPEDANCE


@unit_validator(VOLTAGE)
def divider_voltage_out(thermistor_res: q, basis_res: q) -> q:
    """ Calculates the voltage out of the divider """
    return SUPPLY_VOLTAGE * (basis_res / (thermistor_res + basis_res))


@unit_validator(CURRENT * VOLTAGE)
def self_heating(basis: q, thermistor: q) -> q:
    """ Calculates the self heating of the thermistor """
    return (SUPPLY_VOLTAGE / (basis + thermistor)) ** 2 * thermistor


@unit_validator(CURRENT * VOLTAGE)
def basis_heating(basis: q, thermistor: q) -> q:
    """ Calculates the self heating of the thermistor """
    return (SUPPLY_VOLTAGE / (basis + thermistor)) ** 2 * basis

# @unit_validator(IMPEDANCE)
# def input_impedance(basis: q, thermistor: q) -> q:
#     """ Calculates the input impedance; Assumes no parasitic ind or cap """
#     parallel = (basis * thermistor) / (basis + thermistor) + MUX_RESISTANCE
#     cap_reactance = 1 / (2 * pi * FREQUENCY * (DECODING_CAP + MUX_CAPACITANCE))
#     return (parallel ** 2 + cap_reactance ** 2) ** 0.5


@unit_validator(TIME)
def time_constant(basis: q, thermistor: q) -> q:
    """ Calculates the time constant of the network """
    res = (basis * thermistor) / (basis + thermistor) + MUX_RESISTANCE
    cap = DECODING_CAP + MUX_CAPACITANCE
    return res * cap


def operative_function(basis_res: q) -> tuple:
    """ 
    Returns the voltage range and step_range for a specific basis resistor 
    """
    values = []
    for resistance in THERMISTOR_LOOKUP_TABLE:
        if OPERATING_RANGE[0] > resistance > OPERATING_RANGE[1]:
            values.append(
                divider_voltage_out(resistance, basis_res)
            )

    # Calculates the voltage range and ADC step range
    voltage_range = values[-1]-values[0]
    step_range = (voltage_range) / SUPPLY_VOLTAGE * ADC_RANGE

    return voltage_range, step_range

best_resistor = None
best_ranges = None
#best_impedance = None
best_time_constant = None
basis_heat = None

best_value = float("-inf")

for basis_resistor in standard_resistor_series:
    voltage_range, step_range = operative_function(basis_resistor)
    
    # Calculates the self heating of the thermistor
    max_thm_heating = self_heating(basis_resistor, OPERATING_RANGE[1])
    max_basis_heating = basis_heating(basis_resistor, OPERATING_RANGE[1])
    # Calculates input impedance at min temp; hence max resistance
    # input_imp = input_impedance(basis_resistor, OPERATING_RANGE[0])
    
    # Calculates the time constant at min temp; hence max resistance
    input_time_constant = time_constant(basis_resistor, OPERATING_RANGE[0])
    
    if max_thm_heating > MAXIMUM_SELF_HEAT:
        continue    

    # Calculates that resistor "value"
    value = (
        (step_range / ADC_RANGE) - (max_thm_heating / MAXIMUM_SELF_HEAT) - 
        0.1 * (basis_resistor / max(standard_resistor_series)) - 
        0.5 * (input_time_constant / IDEAL_TIME_CONSTANT)
    )
    
    print(
        f"Resistance: {basis_resistor}, ADC Range: {step_range}, "
        f"Time_constant: {input_time_constant} Basis heating: {max_basis_heating}"
    )
    
    if best_ranges is None:
        best_ranges = [voltage_range, step_range]
        best_resistor = basis_resistor
        # best_impedance = input_imp
        best_time_constant = input_time_constant
        best_value = value
        basis_heat = max_basis_heating
        continue
    
    if value > best_value:
        best_ranges = [voltage_range, step_range]
        best_resistor = basis_resistor
        # best_impedance = input_imp
        best_time_constant = input_time_constant
        best_value = value
        basis_heat = max_basis_heating
        continue
    
print(
    f"Best Resistance: {best_resistor}, Ranges: {best_ranges}, "
    f"best time constant: {best_time_constant} Basis heating: {basis_heat}"
)