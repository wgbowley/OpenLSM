
"""
Filename: low_dropout_regulator.py

Description:
    Computes the power usage by the board and 
    also computes the power loss through the LD117
"""

from picounits import VOLTAGE, TEMPERATURE, CURRENT, POWER, MILLI, MICRO


# LD117 Current Output and thermal limits
input_voltage = 5.0 * VOLTAGE
supply_voltage = 3.3 * VOLTAGE
sot_233_thermal = 110 * TEMPERATURE / POWER

max_quiescent_current = 10 * MILLI * CURRENT
max_output_current = 800 * MILLI * CURRENT
max_output_power = max_output_current * supply_voltage

print(f"Output Power: {max_output_power:.3f}")

# Power Drains
encoder_board = 75.6 * MILLI * POWER

thermistor = 5.496 * MILLI * POWER
thermistor_array = 8 * thermistor

rs422_rs485 = 2.2 * MILLI * CURRENT * supply_voltage
accelerometer = 140 * MICRO * CURRENT * supply_voltage

multiplexer = 16 * MICRO * CURRENT * supply_voltage
stm32 = 30 * MILLI * CURRENT * supply_voltage

power = encoder_board + thermistor_array + rs422_rs485 + accelerometer + multiplexer + stm32
print(f"Total power draw from: {power:.3f}, Ratio: {power/max_output_power*100:.3f} %")

current = power / supply_voltage
power_loss = (input_voltage - supply_voltage) * current + (max_quiescent_current * input_voltage)
total_power = power_loss + power

print(f"Current: {current:.3f}, Power Loss: {power_loss:.3f}, Total Power: {total_power:.3f}")
print(f"Temperature raise: ~{sot_233_thermal * power_loss:.3f}")
