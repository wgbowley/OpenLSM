### 03_boards

#### [00_bridge_driver](00_bridge_driver/readme.md) — Triple Half-Bridge Driver (WIP)

<div align="center">
  <table>
    <tr>
      <td><img src="00_bridge_driver/05_media/mock_up_board_forward.png" alt="Driver Board TOP" style="max-width:375px;"></td>
      <td><img src="00_bridge_driver/05_media/mock_up_board_back.png" alt="Driver Board Bottom" style="max-width:375px;"></td>
    </tr>
    <tr>
      <td><em>Top-side</em></td>
      <td><em>Bottom-Side</em></td>
    </tr>
  </table>
</div>
An isolated triple half-bridge driver with an MCU-side domain of `24 V (DC)` and a power domain of `12-96 V (RMS)`, with current up to `20 A (RMS)`. It supports `step/dir` and `CANBUS` input interfaces and uses `RS-485/RS-422` for communication with external input boards for encoders, Hall-effect sensors, etc. 

---

#### [01_armature_board](01_armature_board/readme.md) — Sensor Board (Manufacturing)

<div align="center">
  <table>
    <tr>
      <td><img src="01_armature_board/05_media/kicad_top_layer.png" alt="Top layer" style="max-width:400px;"></td>
      <td><img src="01_armature_board/05_media/kicad_bottom_layer.png" alt="Bottom layer" style="max-width:400px;"></td>
    </tr>
    <tr>
      <td><em>Top layer — STM32, thermistor array, accelerometer & encoder</em></td>
      <td><em>Bottom layer — base resistors & decoupling capacitors for thermistors</em></td>
    </tr>
  </table>
</div>

A sensor board that connects to the driver over `RS-485/RS-422`. 
It features `8` NTC thermistors switched via an analog multiplexer to produce a thermal profile `T(z, t)`. 
The board also includes a `3-axis SPI accelerometer` and processes data from the encoder board.

---

#### [02_magnetic_encoder](02_magnetic_encoder/readme.md) — Linear Encoder Breakout (Manufacturing)

<div align="center">
  <table>
    <tr>
      <td><img src="02_magnetic_encoder/05_media/kicad_top_layer.png" alt="Top layer" style="max-width:400px;"></td>
      <td><img src="02_magnetic_encoder/05_media/kicad_bottom_layer.png" alt="Bottom layer" style="max-width:400px;"></td>
    </tr>
    <tr>
      <td><em>Top layer — facing the magnetic scale</em></td>
      <td><em>Bottom layer — JST XH 5-pin & mounting side</em></td>
    </tr>
  </table>
</div>

A breakout and mechanical interface for the `AS5311` linear encoder, providing a resolution of `1.95 µm` and an estimated accuracy of `10–20 µm`, depending on mechanical implementation factors. The encoder operates from `0–600 mm/s` and is an incremental `A`, `B`, `Index` type.

---
