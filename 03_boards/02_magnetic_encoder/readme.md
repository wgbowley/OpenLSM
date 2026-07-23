## Overview

OpenLSM uses `closed-loop control` to position the motors armature. This requires a device to measure the armature's position: an encoder. There are many types, magnetic, optical, etc. but OpenLSM uses magnetic linear encoders.

These consist of a linear scale and an encoder head. The encoder measures the magnetic flux density over a `pole pair` (N|S) of the scale to determine its location. This data is output as `A` and `B` binary signals, `90°` out of phase with each other, enabling quadrature decoding `(4× resolution)`. The encoder also outputs a pulse at every `pole pair` transition.

## High-level Topology

The `encoder board`'s main function is to align the encoder mechanically with the magnetic scale. Beyond that, the board sends the `A`, `B`, and `Index` lines to an `STM32` located on the `armature` board.

```
Interface (JST XH 5-pin 2.5mm Pitch Male Header)
Armature Digital Inputs (Digital 3.3 V, STM32)
                    ↓ 
                        ↑
Interface (JST XH 5-pin 2.5mm Pitch Female)
Encoder Board (3.3 V Logic/Domain) 
--------------------------------------------
Magnetic Linear Encoder 
(Power Usage: Typ. 52.8mW - Max. 75.6mW)
(AS5311-ATST-500 - Incremental Mode | A, B, I)
--------------------------------------------
```



## Sampling frequencies

The `AS5311` has a pole pitch of `2 mm` and an A/B pulse ratio of `256`. With quadrature decoding, this ratio is multiplied by `4`, giving `1024` pulses per pole pitch and a step size of `1.96 µm`. It is assumed that to measure the A, B, and Index inputs accurately, the armature board should oversample at `10×` the source frequency.

| Velocity (m/s) | Fundamental (Hz) | Fundamental A/B (Hz) | Sampling (Hz) | Sampling A/B (Hz) |
|:---:|:---:|:---:|:---:|:---:|
| 0.050 | 25.000 | 25.600 k | 250.000 | 256.000 k |
| 0.100 | 50.000 | 51.200 k | 500.000 | 512.000 k |
| 0.150 | 75.000 | 76.800 k | 750.000 | 768.000 k |
| 0.200 | 100.000 | 102.400 k | 1.000 k | 1.024 M |
| 0.250 | 125.000 | 128.000 k | 1.250 k | 1.280 M |
| 0.300 | 150.000 | 153.600 k | 1.500 k | 1.536 M |
| 0.350 | 175.000 | 179.200 k | 1.750 k | 1.792 M |
| 0.400 | 200.000 | 204.800 k | 2.000 k | 2.048 M |
| 0.450 | 225.000 | 230.400 k | 2.250 k | 2.304 M |
| 0.500 | 250.000 | 256.000 k | 2.500 k | 2.560 M |
| 0.550 | 275.000 | 281.600 k | 2.750 k | 2.816 M |
| 0.600 | 300.000 | 307.200 k | 3.000 k | 3.072 M |
| 0.650 | 325.000 | 332.800 k | 3.250 k | 3.328 M |

> [!important]
> This sampling frequency table only covers up to `650 mm/s` as it is the `AS5311`'s maximum linear travel speed.

## Mechanical Considerations

<div align="center">
  <table>
    <tr>
      <td><img src="04_media/top_layer.png" alt="Top layer" style="max-width:400px;"></td>
      <td><img src="04_media/bottom_layer.png" alt="Bottom layer" style="max-width:400px;"></td>
    </tr>
    <tr>
      <td><em>Top layer — Facing the magnetic scale</em></td>
      <td><em>Bottom layer — JST XH 5-pin & mounting side</em></td>
    </tr>
  </table>
</div>

The `linear scale` and `encoder head` should be aligned with each other's centre lines within `~0.5 mm`, and the vertical gap between them should be `~0.3 mm`. The `linear scale` should have a `1 mm` pole length and `2 mm` pole pitch, with a maximum surface flux of `40 mT`. The ambient temperature for a tolerance of `±10 µm` is `-30°C to 70°C`, and the maximum range is `-40°C to 125°C`.

The `encoder board` is `30 mm` in length and `40 mm` in height. The board has a GND and power plane, with the traces embedded into the power plane layer. The board also has M2 bolt holes with a diameter of `~2.10 mm` in a rectangular mounting pattern of `23 mm` and `16.5 mm`. Each M2 bolt hole is directly connected to the GND plane with exposed conductive material around the entire hole.

> [!important]
> The magnetic scale is currently centred on the die centre line. Offsetting the scale relative to the die is not required here, as the `~10 mm` scale width provides sufficient margin for alignment tolerances.

## Documentation

Design notes and implementation decisions are documented in [issue #10](https://github.com/wgbowley/OpenLSM/issues/10).