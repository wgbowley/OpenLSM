### Simulations

#### [00_analytical](./00_analytical/readme.md) — Analytical Tubular Linear Motor Model

<div align="center">
  <img src="../05_media/00_simulation/00_analytical/example.png" alt="Analytical model" style="max-width: 600px">
  <p><em>1D field approximation and FOC showing position (linear) vs force (linear).</em></p>
</div>

This analytical model uses inverse Clarke and Park transforms to compute the phase current based on position, 
then uses 1D field approximations to compute the magnetic co-energy, and finally uses its spatial derivative 
over the z-axis to compute force.

---

#### [00_hybrid](./01_hybrid/readme.md) — Hybrid Tubular Linear Motor Model

*(TBD) — Work in progress*

---
