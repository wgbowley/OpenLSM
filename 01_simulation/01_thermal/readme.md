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