### 06_dropped

This folder contains deliverables that were dropped due to not being on the critical path of `OpenLSM`. However, they are archived here to serve as a starting point for future work, or for others to expand upon.

---

#### [00_bridge_driver](00_bridge_driver) — Triple Half-Bridge Driver

> This item was dropped because a custom bridge driver isn't required, as the `TeensySFOC` works for motor development. This should be continued as a separate repo.

*(Work in progress). The current sensors should be replaced with a magnetic sensor, and the CANBUS interface needs fixing. The schematic needs a design review also.*

A `96 V`, `30 A` triple half-bridge driver with CAN, RS-485, and `step/dir` interfaces, and full isolation between the `STM32`, the half-bridge drivers, and current sensors.

---

#### [01_hybrid](01_hybrid) — Hybrid Tubular Linear Motor Model

> This item was dropped because the `analytical model` works reasonably well as an approximation within appropriate limits.

*(Work in progress). IFEMM works for the FEMM solution and the Biot–Savart kernel works. Just need to combine them and finish the model.*

A hybrid model using `FEMM` to compute the stator field, then using the Biot–Savart law with a parametric slot geometry to compute the armature field. The magnetic co-energy is then calculated assuming uniform magnetic permeability, and finally the spatial derivative over the z-axis is used to compute force.

---