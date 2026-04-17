"""
Filename: thermistor_basis_resistor.py

Description:
    Calculates the basis resistor for the 10kOhm thermistor
    divider circuit into the 74HC4051 and than ADC
"""

from math import pi

from picounits.core import quantities as q, unit_validator
from picounits.constants import (
    IMPEDANCE, VOLTAGE, CURRENT, MILLI, NANO, CAPACITANCE, TIME, PICO
)

from driver.armature_data_board.calculations.thermistor.standard_sizes import (
    standard_resistor_series
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