### Overview

*(TBD) — Work in progress*

---

### Driver & Switches 

Uses the `IRS2104` gate driver for each half bridge and the `SUD90330E` n-channel mosfet for the switches.

*(TBD) — Work in progress*

---

### Power topology

```
Control Side:
INPUT (24v) -> LM2596S-12 (12v) LM2596S-5 (5v) -> LD1117V33 (3.3v) 

------------------------------------------------------------------

Isolated Side:
Input (up to 96V), Isolated 12V & Isolated 3.3V
```

---

### Documentation
Reference: [Design and build triple bridge driver based on SimpleFOC shield](https://github.com/wgbowley/OpenLSM/issues/12)

---