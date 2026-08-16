## Overview

<div align="center">
  <table>
    <tr>
      <td><img src="04_media/top_layer.png" alt="Top layer" style="max-width:400px;"></td>
      <td><img src="04_media/bottom_layer.png" alt="Bottom layer" style="max-width:400px;"></td>
    </tr>
    <tr>
      <td><em>Top layer — STM32, Thermistor Array, Accelerometer & Encoder</em></td>
      <td><em>Bottom layer — Basis resistors & decoupling caps for Thermistors</em></td>
    </tr>
  </table>
</div>


OpenLSM uses `closed-loop control` to position the motor's armature. This requires a device to measure the armature's position: [an encoder board](../02_magnetic_encoder). OpenLSM also requires acceleration and thermal data for validating the motor's transient behaviour during operation. To achieve this, an `SPI` accelerometer is used to collect `3-axis` acceleration data, and an array of `8` thermistors across the motor's `z-axis` is used to collect temperature data `T(z, t)`.

The array consists of `8` NTC thermistors and an analog multiplexer, which feeds into the ADC pin on the peripheral controller. The data collected from these sensors is then processed into relevant quantities and sent over the `RS-485/RS-422` link to the motor controller board.

## High-level topology

The `armature board`'s main function is to collect the data from the integrated accelerometer, integrated thermal array and encoder board, then process the data into quantities and send it to the motor controller board. Given the locality to the armature slots and their tendency to produce electromagnetic noise, an `RS-422/RS-485` differential link was used.

```
Interface (JST XH 4-pin 2.5mm Pitch Male Header)
Motor Controller Board (Digital 3.3 V, 5V STM32)
↓
↑
Interface (JST XH 4-pin 2.5mm Pitch Female)
Armature Board (3.3 V Logic/domain)
-------------------------------------------------
Low dropout regulator (5V - 3.3V) (LD1117) (800mA Max Current) (~20°C Expected)
↓
RS-485/RS-422 Link (MAX3485) (Differential Pair Standard) or (Rx, Tx Bypass)
↓ ↑
Peripheral Controller (170 MHz) (STM32G431K8Tx)
↓ ↑
Thermistor Array (800 Hz), SPI Accelerometer, Linear Encoder Interconnect (350 kHz)
-------------------------------------------------
```


## Mechanical Considerations

<div align="center">
  <table>
    <tr>
      <td><img src="04_media/top_layer.png" alt="Top layer" style="max-width:400px;"></td>
      <td><img src="04_media/bottom_layer.png" alt="Bottom layer" style="max-width:400px;"></td>
    </tr>
    <tr>
      <td><em>Top layer — STM32, Thermistor Array, Accelerometer & Encoder</em></td>
      <td><em>Bottom layer — Basis resistors & decoupling caps for Thermistors</em></td>
    </tr>
  </table>
</div>


The `armature board` is intended to be located either tangentially to the armature or below the armature in parallel with the linear rail carriage face. The board is `60 mm` in length and `40 mm` in width. The board has a GND and power plane, with the traces embedded into the power plane layer. The board also has M2 bolt holes with a diameter of `~2.10 mm` in a rectangular mounting pattern of `54 mm` and `34 mm`. Each M2 bolt hole is directly connected to the GND plane with exposed copper around the entire hole.

## Software considerations

*(TBD) — Work in progress*

> [!important]
> Programming is done via a `2×3 pin`, `2.54 mm` vertical male connector on the armature board, located near the `STM32G431K8Tx`.

## Documentation

Design notes and implementation decisions are documented in [issue #11](https://github.com/wgbowley/OpenLSM/issues/11).