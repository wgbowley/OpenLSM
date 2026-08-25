### 00_analytical

> [!WARNING]
>
> Do not assume the resulting design variables are suitable for fabrication or real-world use without further analysis.

This analytical model uses 1D field approximations and inverse- clark, park transformers to move the motor armature 
along the stator to get a force vs position curve.

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