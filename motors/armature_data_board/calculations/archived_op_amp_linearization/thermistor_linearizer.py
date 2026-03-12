"""
Filename: thermistor_linearizer.py

Description:
    Calculates the resistors needed to produce a linear output voltage
    between 0 and 3.3 volts for the pi pico ADC line.
"""

from math import exp

from picounits.core import unit_validator, Quantity as q
from picounits import (
    IMPEDANCE, VOLTAGE, TEMPERATURE, DIMENSIONLESS, KILO, MILLI
)

""" Thermistor Steinhart-Hart Inverse Coefficients """
A1 = -14.6337
B1 = 4791.84
C1 = -115334
D1 = -3.7305e+06

""" Circuit Variables """
SUPPLY = 3.3 * VOLTAGE
V_OUT_MAX = 3.2 * VOLTAGE
V_OUT_MIN = 100 * MILLI * VOLTAGE

MAX_TEMPERATURE = 125 * TEMPERATURE
MIN_TEMPERATURE = 25 * TEMPERATURE

NTC100AT25C = 100 * KILO * IMPEDANCE

""" Standard value for the feedback resistor """
FEEDBACK = 1.5 * KILO * IMPEDANCE


@unit_validator(IMPEDANCE)
def calculate_ntc_resistance(temp_celsius: q) -> q:
    """ Calculates the NTC thermistor resistance using Steinhart-Hart """
    tk = temp_celsius.value + 273.15
    exponent_value = (A1 + (B1 / tk) + (C1 / tk ** 2) + (D1 / tk ** 3))

    return NTC100AT25C * exp(exponent_value)


""" Calculated values for resistance at max temp and min temp """
R_NTC_MAX = calculate_ntc_resistance(MAX_TEMPERATURE)
R_NTC_MIN = calculate_ntc_resistance(MIN_TEMPERATURE)


@unit_validator(IMPEDANCE)
def ntc_low_side_resistor() -> q:
    """ Calculates the low side resistor of the divider """
    return (R_NTC_MIN * R_NTC_MAX) ** 0.5


@unit_validator(VOLTAGE)
def voltage_range(low_side_ntc) -> tuple[q]:
    """ Calculates the input voltage range, minium and maximum """
    def calculation(temp_res):
        """ Calculation for voltage depends on different temp res"""
        return SUPPLY * (low_side_ntc / (temp_res + low_side_ntc))

    return calculation(R_NTC_MIN), calculation(R_NTC_MAX)


@unit_validator(DIMENSIONLESS)
def ideal_operational_gain(v_in_min: q, v_in_max: q) -> q:
    """ Calculates the ideal operational gain required for output swing """
    return (V_OUT_MAX - V_OUT_MIN) / (v_in_max - v_in_min)


@unit_validator(IMPEDANCE)
def parallel_resistance_ref_divider(ideal_gain: q) -> q:
    """ Calculates the parallel resistance of the two elements """
    return FEEDBACK / (ideal_gain - 1)


@unit_validator(IMPEDANCE)
def ref_divider_resistance(
    v_in_max: q, ideal_gain: q, parallel: q
) -> tuple[q]:
    """ Calculates the resistance of the high_side and low_side resistors """
    high_side = FEEDBACK * SUPPLY / (v_in_max * ideal_gain - V_OUT_MAX)
    low_side = parallel * high_side / (high_side - parallel)
    return high_side, low_side


@unit_validator(DIMENSIONLESS)
def actual_gain(parallel: q) -> q:
    """ Calculates the actual op-amp gain of the system"""
    return (parallel + FEEDBACK) / parallel


""" Calculations and outputs """
ntc_low = ntc_low_side_resistor()
v_in_min, v_in_max = voltage_range(ntc_low)
ideal_gain = ideal_operational_gain(v_in_min, v_in_max)
parallel = parallel_resistance_ref_divider(ideal_gain)
ref_high, ref_low = ref_divider_resistance(v_in_max, ideal_gain, parallel)
gain = actual_gain(parallel)

print("============================")
print(f"Resistance at {MIN_TEMPERATURE.value} c is {R_NTC_MIN:.3f}")
print(f"Resistance at {MAX_TEMPERATURE.value} c is {R_NTC_MAX:.3f}")
print("============================")

# Print out names based partly on naming in the schematic present in (SBOA323A)
print(f"Low side NTC divider (R1) = {ntc_low:.3f}")
print(f"Low side reference divider (R2) = {ref_low:.3f}")
print(f"High side reference divider (R3) = {ref_high:.3f}")
print("============================")
print(f"Op Amp feedback (R4) = {FEEDBACK:.3f}")
print(f"Op Amp gain = {gain:.3f}")
print("============================")
