### 00_analytical

This analytical model uses inverse Clarke and Park transforms to compute the phase current based on position, then uses 1D field 
approximations to compute the magnetic co-energy, and finally uses its spatial derivative over the z-axis to compute force.

#### High Level Topology

```
Z-axis Position (z)
    ↓
Electrical Angle (Radians)
    ↓
Phase Currents (I_a, I_b, I_c)
    ↓
1D Magnetostatic Fields (Armature, Stator)
    ↓
Integration of Fields for Co-Energy
    ↓
Spatial Derivative of Co-Energy Along Z-axis
    ↓
Force on Armature at That Position
    ↺ (Moves to next position)
```

> Do not assume the resulting design variables are suitable for fabrication or real-world use without further analysis.

---

### Mathematical Implementation

#### Magnetic Force Derivation

It is assumed that the field dynamics within a tubular linear motor can be compressed into 1D representations. With these fields, the force can be calculated using virtual energy methods:

$$ F = -\frac{dU}{dz} \approx -\frac{U(z+\Delta z) - U(z-\Delta z)}{2\Delta z} $$

The magnetic energy distributed through the motor is calculated from the integral of the $B$ and $H$ fields:

$$ U = \frac{1}{2} \int_{\text{vol}} (H_{\text{stator}} + H_{\text{armature}}) \cdot (B_{\text{stator}} + B_{\text{armature}}) \, dv $$

Due to the relationship between the fields being $B = \mu H$ and $\mu_r = 1$ for ironless motors, the integral can be simplified:

$$ U = \frac{\mu_0}{2} \int_{\text{vol}} (H_{\text{stator}} + H_{\text{armature}})^2 \, dv $$

$$ U = \frac{\mu_0}{2} \int_{\text{vol}} (H_{\text{stator}}^2 + 2 H_{\text{armature}} H_{\text{stator}} + H_{\text{armature}}^2) \, dv $$

Because the spatial derivative is being calculated, both $H_{\text{stator}}^2$ and $H_{\text{armature}}^2$ can be removed as they are constants:

$$ U = \mu_0 \int_{\text{vol}} (H_{\text{armature}} H_{\text{stator}}) \, dv $$

The volume integral can then be simplified for the 1D approximation using $A_{\text{rad}}$ between the slot mean and pole outer diameter:

$$ U = A_{\text{rad}} \cdot \mu_0 \int_{\Omega} H_{\text{armature}}(z - z_{o}) \cdot H_{\text{stator}}(z) \, dz $$

Where $z_o$ is the armature offset with respect to the stator, and $\Omega$ is the integration domain over the z-axis.


#### 1D Stator Approximation

The stator field is represented using an analytical approximation based on hyperbolic-secant field distributions, which were derived from equivalent finite element solutions:

$$
H_{\text{dipole}}(z, m) = h_{\text{peak}} \left[ \mathrm{sech}\left( \frac{n\left((z-m) +\ell_{\text{dipole}}/2 \right)}{\ell_{\text{dipole}}} \right)^2 -\mathrm{sech}\left(\frac{n\left((z-m) -\ell_{\text{dipole}}/2 \right)}{\ell_{\text{dipole}}} \right)^2 \right]
$$

The parameter $n$ determines the decay rate of the dipole field along the z-axis, and $m$ is the translation along the z-axis. This formula is then used as a dipole primitive to produce an `N-S-S-N` stator field:

$$H_{\text{stator}}(z) = \sum_{i = 0}^{N_{\text{dipoles}}} (-1)^i \cdot H_{\text{dipole}}(z, \, i \cdot L_{\text{dipole}})$$

where $N_{\text{dipoles}}$ is the number of dipoles within the stator and $L_{\text{dipole}}$ is the length of an individual pole.


#### 1D Armature Approximation

The armature field is represented using the finite-length Biot–Savart formulation:

$$
H_{\text{pole}}(z) = \frac{N I}{2\ell} \int_{-\ell_{\text{pole}}/2}^{\ell_{\text{pole}}/2} 
\frac{R^2}{\left[(z - z')^2 + R^2\right]^{3/2}} \, dz'
$$

where $N$ is the number of turns within a slot and $I$ is the current moving through that slot's phase. This field is then used as a pole primitive to produce the armature fields via superposition:

$$H_{\text{armature}}(z) = \sum_{i=0}^{N_{\text{slots}}} (-1)^{i+1} \cdot h_{\text{coil}}\left(z, \, z_i, \, I_{\text{ph}(i)}, \, N_{\text{turns}}\right)$$

where $N_{\text{slots}}$ is the number of slots and $(-1)^{i+1}$ encodes the slot polarity within the armature.

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
Solves the solution (virtual work methods)
        ↓
============================================
Re-attaches units at validated boundaries
============================================
        ↓
Generates matplotlib graph (force vs position)
--------------------------------------------
        ↓
Output (matplotlib + terminal)
--------------------------------------------
matplotlib GUI + Terminal printout
--------------------------------------------
```

#### Numerical Flow

The model samples the motor over a displacement range. At each position, the electrical angle is calculated, phase currents are updated, the armature and stator fields are evaluated, and their interaction is integrated along the z-axis.

The force is obtained using a central finite-difference approximation of co-energy and
the integration domain is terminated when the field interaction remains below a specified
threshold for a fixed number of samples.

```py
while True:
    # Computes field strength of pole & coil at z
    h_armature = self._armature_field(z, translate)
    h_stator = self._stator_field(z)

    # Computes the co-energy from that interaction
    interaction = h_armature * h_stator
    interactions.append(interaction)

    if abs(interaction) <= epsilon:
        below_epsilon += 1
    else:
        below_epsilon = 0

    # Breaks loop if below epsilon for more than window iterations
    if below_epsilon >= window:
        break

    # Moves along the z-axis
    z += dz
```

---

### Model Validation

*Work in progress*

The analytical model is intended as a reduced-order model rather than a replacement for finite-element analysis (FEA). It is quite quick, though exact speeds depend on parameters. It is useful for understanding the mechanics of linear motors without the large execution times of FEA.

---

### Documentation

Reference: [Build a reduced order model for tubular linear synchronous motors](https://github.com/wgbowley/OpenLSM/issues/7)

---