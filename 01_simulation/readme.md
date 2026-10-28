### 01_simulations

#### [00_magnetic](./00_magnetic/readme.md) - Magnetic Model

This analytical model uses control currents `(id, iq)` and then a three-phase transform to compute the phase current based on position, 
then uses 1D field approximations to compute the magnetic co-energy, and finally uses its spatial derivative over the z-axis to compute force. 
This is then repeated across the stator to produce a force vs position graph.


<div align="center">
  <img src="../05_media/00_simulation/00_magnetic/example.png" alt="Analytical model" width="700">
  <p><em>1D field approximation and FOC showing position (linear) vs force (linear).</em></p>
</div>

---

#### [01_thermal](./01_thermal/readme.md) - Thermal Model

This analytical model uses heat conduction to model the temperature drop across the potting material and heat sink base. It then uses the surface temperature
of the exposed heat sink to calculate the conductive, convective and radiate heat dissipation of a specific geometry. It then repeats this across a grid
of fin spacing and radial fin sizes to produce a thermal performance graph.

<div align="center">
  <img src="../05_media/00_simulation/01_thermal/example.png" alt="Thermal Model" width="700">
  <p><em>Grid search over fin axial spacing and outer diameter.</em></p>
</div>

---