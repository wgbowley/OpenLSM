### 01_thermal

This analytical model is intended as a reduced-order model rather than a replacement for computational thermal fluid dynamics. It is rather quick, however, at sampling the solution space. It is useful for understanding the parameters of the annulus heat sink without long execution times.

#### High Level Topology

> Do not assume the resulting design variables are suitable for fabrication or real-world use without further analysis.

```
Temperature delta across potting and heat-sink base
    ↓
Rayleigh Number (Flow fluid Regime)
    ↓
Elenbaas Correlation for Natural Convection
    ↓
Modified Bessel Function (r1, r2 fixed temperatures)
    ↓
Geometric Fin Convection Efficiency
    ↓
Total Fin and Barrel Area based on Number of Fins
    ↓
Radiative Dissipation via Stefan-Boltzmann
    ↓
Convection via Fins and Barrel Area Together
    ↓
Total Dissipation via Radiative and Convection Dissipation
```

> *(Note). Temperature delta is only calculated during initialization. It is assumed fixed.* 

---

### Mathematical Implementation

*(TBD) - Work in progress*

#### Conduction Rate (Fourier's Law)

Fourier's law is used to compute the temperature drop across the potting material around the coils, and then the temperature drop across the annulus heat-sink base. It is assumed that thermal energy is conserved and migrating outwards.

$$ Q = -kA\nabla T \rightarrow -kA\frac{dT}{dr} \Rightarrow dT = -\frac{Q}{kA} \, dr $$

The differential equation can then be solved assuming the area is $A = 2\pi r L$:

$$ \int_{T_1}^{T_2} dT = -\int_{r_1}^{r_2} \frac{Q}{k(2\pi r L)} \, dr \Rightarrow T_2 - T_1 = -\frac{Q}{2\pi k L} \ln\left(\frac{r_2}{r_1}\right) $$

The final formula for the temperature drop, or more precisely the temperature delta across the cylindrical layer:

$$ \Delta T = -\frac{Q}{2\pi k L} \ln\left(\frac{r_2}{r_1}\right) $$

---

### Computational Implementation

#### High-Level Topology

```
UnitValues
--------------------------------------------
parameters.uiv ← derived.ut
--------------------------------------------
        ↓
Python +3.10
--------------------------------------------
Solver Class
        ↓
Extracts Qualities from parameters
        ↓
============================================
Validates & removes units from parameters
============================================
        ↓
Solves the performance for spacing/thickness
        ↓
--------------------------------------------
        ↓
Generates matplotlib graph (dissipation (W))
--------------------------------------------
        ↓
Output (matplotlib + terminal)
--------------------------------------------
matplotlib GUI
--------------------------------------------
```

---

### Documentation

Reference: [`Refine the annulus heat sink model for tubular linear motors into a 01_simulation model.`
](https://github.com/wgbowley/OpenLSM/issues/60)

---