### Armature Board — Bill of Materials (BOM)

The procurement file [`04_procurement/digikey.csv`](./04_procurement/digikey.csv) is for `2x` 
boards and does not include the JST socket or the `2x03` male header.

---

| Status | Reference | Qty | Value | Footprint | Notes / Datasheet |
|--------|-----------|-----|-------|-----------|-------------------|
| [x] | BOOT1, RESET1 | 2 | SW_Push | Button_Switch_SMD:SW_Push_1P1T_NO_CK_KMR2 | |
| [x] | C1, C6 | 2 | 1µF | Capacitor_SMD:C_0603_1608Metric | 10 V rated |
| [x] | C3, C5, C9, C13, C14, C15 | 6 | 100nF | Capacitor_SMD:C_0603_1608Metric | 16 V rated |
| [x] | C4, C7 | 2 | 10µF | Capacitor_SMD:C_1206_3216Metric | 10 V rated |
| [x] | C8 | 1 | 22µF | Capacitor_SMD:C_1206_3216Metric | 25 V rated |
| [x] | C10, C16, C17, C18, C19, C20, C21, C22, C23 | 9 | 10nF | Capacitor_SMD:C_0603_1608Metric | 16 V rated |
| [x] | C11, C12 | 2 | 30pF | Capacitor_SMD:C_0603_1608Metric | 50 V rated |
| [x] | D1 | 1 | LED | LED_SMD:LED_0603_1608Metric_Pad1.05x0.95mm_HandSolder | |
| [-] | ENC1 | 1 | 01x05 | Connector_JST:JST_PH_B5B-PH-K_1x05_P2.00mm_Vertical |  Any standard JST socket set. |
| [x] | FB1 | 1 | BLM18AG601SN1D | BLM18AG601SN1D:BEADC1608X95N | |
| [-] | INT1 | 1 | 01x04 | Connector_JST:JST_PH_B4B-PH-K_1x04_P2.00mm_Vertical |  Any standard JST socket set. |
| [-] | NTC1-NTC8 | 8 | 01x02 | Connector_JST:JST_PH_B2B-PH-K_1x02_P2.00mm_Vertical |  Any standard JST socket set. |
| [-] | PRO1 | 1 | 01x06 | Connector_PinHeader_2.54mm:PinHeader_2x03_P2.54mm_Vertical | Any standard 2.54 set. |
| [x] | R1, R2 | 2 | 10Ω | Resistor_SMD:R_0603_1608Metric | |
| [x] | R3 | 1 | 120Ω | Resistor_SMD:R_0603_1608Metric | |
| [x] | R5, R6 | 2 | 10kΩ | Resistor_SMD:R_0603_1608Metric | |
| [x] | R7, R16, R17 | 3 | 220Ω | Resistor_SMD:R_0603_1608Metric | |
| [x] | R8-R15 | 8 | 1.8kΩ | Resistor_SMD:R_0603_1608Metric | +/- 0.1% |
| [x] | U1 | 1 | STM32G431K8Tx | Package_QFP:LQFP-32_7x7mm_P0.8mm | Instead: STM32G431K6T6 |
| [x] | U2 | 1 | LD1117S33 | Package_TO_SOT_SMD:SOT-223-3_TabPin2 | |
| [x] | U3 | 1 | MAX3485 | Package_SO:SOIC-8_3.9x4.9mm_P1.27mm | |
| [x] | U4 | 1 | ADXL343 | Package_LGA:LGA-14_3x5mm_P0.8mm_LayoutBorder1x6y | |
| [x] | U5 | 1 | 74HC4051 | 74HC4051:SSOP-16 | |
| [x] | Y1 | 1 | ECS-80-20-18-TR | ECS-80-18-23B-JTN-TR:XTAL_ECS-80-20-18-TR | |

> *(Note). `[ ]` Not Ordered. `[-]` Not Required. `[x]` Ordered.*

---