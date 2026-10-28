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

$$ \beta = \frac{1}{\Delta T_{fin} + T_{env}} \Rightarrow Ra = \frac{g \beta \Delta T_{fin} \ell^3}{\nu \alpha} $$

Where $\nu$ is kinematic viscosity and $\alpha$ is thermal diffusivity. The Rayleigh-Elenbaas number and correlation were then computed:

$$ Ra^* = Ra \times \frac{\ell}{t} \Rightarrow Nu = \frac{1}{24} Ra^* \left(1 - e^{-\frac{35}{Ra^*}}\right)^{\frac{3}{4}} $$

Where $\ell$ is the axial spacing between fins and $t$ is the radial thickness of the fin. From the Nusselt number, the heat transfer coefficient could then be computed:

$$ h = \frac{Nu \cdot k}{L} $$

Where $L$ is the axial characteristic length of the heat-sink.

#### Fin Efficiency

The efficiency of the fin was computed using the governing heat transfer equation for an annular fin:

$$ z^2 \frac{d^2\theta}{dz^2} + z \frac{d\theta}{dz} - z^2 \theta = 0 $$

It can be observed that this specific differential equation matches the modified Bessel equation with $\alpha = 0$:

$$ z^2 \frac{d^2\theta}{dz^2} + z \frac{d\theta}{dz} - (z^2 - \alpha^2)\theta = 0 $$
$$ z^2 \frac{d^2\theta}{dz^2} + z \frac{d\theta}{dz} - z^2\theta = 0 $$

The parameter $\gamma$ can then be computed for the ratio of conduction to convection:

$$ \gamma = \sqrt{\frac{h}{k w}} $$

where $h$ is the heat transfer coefficient of the environment and $k$ is the thermal conductivity of the fin material. The efficiency parameter can then be calculated using the first and second kind solutions of the Bessel function with fixed boundaries:

$$ u_0 = \gamma \, r_{1}, \quad u_1 = \gamma \, r_{2} $$

$$ t_1 = K_1(u_0) \, I_1(u_1) - I_1(u_0) \, K_1(u_1), \quad t_2 = K_0(u_0) \, I_1(u_1) + I_0(u_0) \, K_1(u_1) $$

$$ \eta = \frac{2 r_1}{\gamma (r_2^2 - r_1^2)} \frac{t_1}{t_2} $$

Where $r_2$ is the outer radius of the fin and $r_1$ is the inner radius of the fin.

#### Convection & Radiate Heat Dissipation

The convection heat dissipation was calculated by reducing the problem to a single fin:

$$ A_f = 4\pi \left( r_2^2 - r_1^2 \right) $$

$$ Q_f = \eta \, h \, A_f \, \Delta T_f $$

The number of fins is then computed using `mod` and `floor` to estimate the count, and ultimately $Q_f$ is multiplied by that estimate. The convection dissipation from the "barrel" of the heat-sink is then computed:

$$ A_b = 2\pi r_1 \ell_{l-f} $$

$$ Q_b = h A_b \Delta T_f $$

The radiative heat dissipation was then computed using the Stefan–Boltzmann law with the total area $A_t = A_b + A_f$:

$$ Q_r = \varepsilon \sigma A_t \left( (T_{env} + \Delta T_f)^4 - T_{env}^4 \right) $$

The total heat dissipation is then simply the sum:

$$ Q_t = Q_f + Q_b + Q_r $$

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