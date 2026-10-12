### Encoder Board — Bill of Materials (BOM)

The procurement file [`04_procurement/digikey.csv`](./04_procurement/digikey.csv) is 
for `2x` boards and does not include the JST socket or the R1 `0 Ω` resistor.

---

| Status | Reference | Qty | Value | Footprint | Notes / Datasheet |
|--------|-----------|-----|-------|-----------|-------------------|
| [x] | C1, C3 | 2 | 100 nF | Capacitor_SMD:C_0402_1005Metric | 10V Rated |
| [x] | C2 | 1 | 10 µF | Capacitor_SMD:C_1206_3216Metric | 10V Rated |
| [-] | J1 | 1 | 01x05 | Connector_JST:JST_XH_B5B-XH-A_1x05_P2.50mm_Vertical | Any standard JST socket set. |
| [-] | R1 | 1 | 0 Ω | Resistor_SMD:R_0402_1005Metric | Can be ignored — Just solder it. |
| [x] | U1 | 1 | AS5311-ATST-500 | footprints:TSSOP20_AS5311-ATST-500_AMS | |

> *(Note). `[ ]` Not Ordered. `[-]` Not Required. `[x]` Ordered.*

---