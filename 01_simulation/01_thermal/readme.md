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

#### Conduction Rate (Fourier's Law)

Fourier's law is used to compute the temperature drop across the potting material around the coils, and then the temperature drop across the annulus heat-sink base. It is assumed that thermal energy is conserved and migrating outwards.

$$ Q = -kA\nabla T \rightarrow -kA\frac{dT}{dr} \Rightarrow dT = -\frac{Q}{kA} \, dr $$

The differential equation can then be solved assuming the area is $A = 2\pi r L$:

$$ \int_{T_1}^{T_2} dT = -\int_{r_1}^{r_2} \frac{Q}{k(2\pi r L)} \, dr \Rightarrow T_2 - T_1 = -\frac{Q}{2\pi k L} \ln\left(\frac{r_2}{r_1}\right) $$

The final formula for the temperature drop, or more precisely the temperature delta across the cylindrical layer:

$$ \Delta T = -\frac{Q}{2\pi k L} \ln\left(\frac{r_2}{r_1}\right) $$

#### Natural Convection & Heat transfer

The heat transfer coefficient of the annulus heat-sink is computed via the Rayleigh-Elenbaas correlation for natural convection. The Rayleigh number was computed assuming constant pressure during thermal expansion:

$$ \beta = \frac{1}{T_{\text{fin}} + T_{\text{env}}} \Rightarrow Ra = \frac{g \beta T_{\text{fin}} \ell^3}{\nu \alpha} $$

Where $\nu$ is kinematic viscosity and $\alpha$ is thermal diffusivity. The Rayleigh-Elenbaas number and correlation were then computed:

$$ Ra^* = Ra \times \frac{\ell}{t} \Rightarrow Nu = \frac{1}{24} Ra^* \left(1 - e^{-\frac{35}{Ra^*}}\right)^{\frac{3}{4}} $$

Where $\ell$ is the axial spacing between fins and $t$ is the radial thickness of the fin. From the Nusselt number, the heat transfer coefficient could then be computed:

$$ h = \frac{Nu \cdot k}{L} $$

Where $L$ is the axial characteristic length of the heat-sink.

#### Fin Efficiency

*(TBD) - Work in progress*

#### Convection & Radiate Heat Dissipation

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