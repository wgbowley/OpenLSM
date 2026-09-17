#### Overview

<div align="center">
  <table>
    <tr>
      <td><img src="05_media/kicad_top_layer.png" alt="Top layer" style="max-width:400px;"></td>
      <td><img src="05_media/kicad_bottom_layer.png" alt="Bottom layer" style="max-width:400px;"></td>
    </tr>
    <tr>
      <td><em>Top layer — Teensy & SimpleFOC shield</em></td>
      <td><em>Bottom layer — Supporting electronics</em></td>
    </tr>
  </table>
</div>

The Teensy SFOC board is a breakout board for the Teensy 4.1 and SimpleFOC Arduino shield, with `RS-485/RS-422` support and a secondary encoder input if `RS-485` to the armature board cannot be used. It also has a `step/dir` input for use with standard 3D printer main-boards.

> *(Note). This is a prototyping/dev board. It is not meant to be used long-term.*

---

### Documentation

Design notes and implementation decisions are documented in [issue #42](https://github.com/wgbowley/OpenLSM/issues/42).

---