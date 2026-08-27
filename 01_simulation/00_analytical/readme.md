### 00_analytical

This analytical model uses inverse Clarke and Park transforms to compute the phase current based on position, then uses 1D field 
approximations to compute the magnetic co-energy, and finally uses its spatial derivative over the z-axis to compute force.

> Do not assume the resulting design variables are suitable for fabrication or real-world use without further analysis.

---

### Mathematical Implementation

Force calculation based on virtual energy methods:

$$ F = -\frac{dU}{dz} \approx -\frac{U(z+\Delta z) - U(z-\Delta z)}{2\Delta z} $$

is calculated from the integral of the energy distributed through the motor:

$$ U = \frac{1}{2} \int_{\text{vol}} (H \cdot B) \, dv $$

---

### Computational Implementation

*(TBD) — Work in progress*

---

### Model Validation

*(TBD) — Work in progress*

---

### Documentation

Reference: [Build a reduced order model for tubular linear synchronous motors](https://github.com/wgbowley/OpenLSM/issues/7)

---