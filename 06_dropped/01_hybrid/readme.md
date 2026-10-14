### 01_hybrid

This hybrid model uses the inverse Clarke and Park transforms to compute the phase current based on position, then uses the Biot–Savart law to compute a single slot sheet, and then duplicates and shifts it to form the slot kernel. That kernel and the `FEMM` solution are used to calculate the magnetic co-energy, and finally the spatial derivative over the z-axis is used to compute force.

#### High Level Topology

*(TBD) - Work in progress*

```

```

> Do not assume the resulting design variables are suitable for fabrication or real-world use without further analysis.

---

### Mathematical Implementation
*(TBD) - Work in progress*

---

### Computational Implementation
*(TBD) - Work in progress*

---

### Model Validation
*(TBD) - Work in progress*

---

### Documentation

Reference: [`Build a reduced order model for tubular linear synchronous motors`](https://github.com/wgbowley/OpenLSM/issues/7)<br>
Reduced Order Modelling tool: [`ifemm`](https://github.com/wgbowley/iFEMM)

---