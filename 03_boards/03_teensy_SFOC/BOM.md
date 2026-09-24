### TeensySFOC — Bill of Materials (BOM)

*(TBD) — Work in progress*

---

| Status | Reference | Qty | Value | Footprint | Notes / Datasheet |
|--------|-----------|-----|-------|-----------|-------------------|
| [ ]   | C1, C2    | 2   | 1uf   | Capacitor_SMD:C_1206_3216Metric | Decoupling |
| [ ]   | C3, C4    | 2   | 100nf | Capacitor_SMD:C_1206_3216Metric | Decoupling |
| [ ]   | D1, D2    | 2   | 1N4148 | Diode_SMD:D_SOD-123 | |
| [ ]   | ENC1      | 1   | 01x05 | Connector_JST:JST_XH_B5B-XH-A_1x05_P2.50mm_Vertical | Encoder input |
| [ ]   | INT1      | 1   | 01x04 | Connector_JST:JST_XH_B4B-XH-A_1x04_P2.50mm_Vertical | Interface |
| [ ]   | PW1, SD1  | 2   | 01x02 | Connector_JST:JST_XH_B2B-XH-A_1x02_P2.50mm_Vertical | Power / Shutdown |
| [ ]   | R1, R2    | 2   | 10    | Resistor_SMD:R_1206_3216Metric | |
| [ ]   | R3        | 1   | 120   | Resistor_SMD:R_1206_3216Metric | |
| [ ]   | R4, R5    | 2   | 1k    | Resistor_SMD:R_1206_3216Metric | |
| [ ]   | SFOC1     | 1   | SimpleFOC Shield | Module:Arduino_UNO_R3 | [SimpleFOC]() |
| [ ]   | U1        | 1   | Teensy4.1 | teensy_library:Teensy41 | Main MCU |
| [ ]   | U2        | 1   | MAX3485 | Package_SO:SOIC-8_3.9x4.9mm_P1.27mm | [Datasheet](https://datasheets.maximintegrated.com/en/ds/MAX3483-MAX3491.pdf) — RS-485 transceiver |
| [ ]   | U3        | 1   | LD1117V33 | Package_TO_SOT_SMD:SOT-223-3_TabPin2 | [Datasheet](https://www.st.com/resource/en/datasheet/ld1117.pdf) — 3.3V LDO regulator |

> *(Note).* `[ ]` Not Ordered. `[-]` Not Required. `[x]` Ordered.

---