<!--
Color Palette:
FFFFFF - pure white 
FF8F0E - bold, warm orange with a strong golden-yellow undertone

I'm going to do a little experiment here and not update the logo
until someone raises an issue about light-mode readability.
- William Bowley, 2026-07-05

P.S: Thanks for downloading the OpenLSM repository `▽`ʃ♡
-->

<p align="center">
  <img src="05_media/01_logos/logo.png" alt="OpenLSM" style="max-width:600px;">
  <br>
  <em>
    Low Cost Linear Synchronous Motors – Designed & built by 
    <a href="https://github.com/wgbowley">William Bowley</a> (primary) &amp; 
    <a href="https://github.com/LawsonDG">Lawson Gallup</a>
  </em>
</p>

## Overview
![Status](https://img.shields.io/badge/Status-WIP-FFFFFF?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-FFFFFF?style=flat-square&color=white)
![Focus](https://img.shields.io/badge/Focus-Simulation-FF8F0E?style=flat-square)
![Domain](https://img.shields.io/badge/Domain-Hardware-FF8F0E?style=flat-square&color=FF8F0E)

OpenLSM is an experimental project with the objective of designing low-cost permanent magnet linear motors for Cartesian motion systems such as pick-and-place machines or CNC machines. The project will fulfill this goal by using readily available materials and tooling, combined with reduced-order and finite element models.

### Objectives

> [!IMPORTANT]
> - Design low-cost permanent magnet linear motors for Cartesian motion systems.
> - Develop computational models that are validated against experimental measurements.
> - Demonstrate continuous thermal steady-state operation under defined operating conditions.
>
> *More specific targets `(Force density, Efficiency, Cost, Etc.)` to be defined as the domain becomes clearer.*

## Methodology

The methodology for openLSM is to design a reduced-order or finite element model to compute the expected results for a specific motor topology. The resulting model is used to tune design parameters either manually or via algorithmic optimization. The parameters are then used to produce a CAD model, which is fabricated and experimentally tested. Discrepancies between predicted and observed performance are used to update the model for subsequent iterations.

```
Conceptual Design ↔ Reduced-order / FEA Model
      ↓
Design Parameters ↔ Detailed Design (CAD)
      ↓
Fabrication & Testing
      ↓
Post-Analysis & Model Validation
      ↺
```

<div align="center">
  <em>
  Figure 1: OpenLSM design methodology. Computational models are iteratively refined through experimental validation.</em>
</div>

## Prototype Alpha

An ironless flat linear motor with a polylactic acid (PLA) armature featuring `6` slots, hand wound using `0.2 mm` diameter enameled copper wire with `2` slots in-series per phase. The stator, similar to the armature, was printed in PLA and had `4` pole pairs per armature length and `10` pole pairs total. The motor produced measurable force, though the magnitude was not quantified before the PLA coil forms deformed due to thermal stress.

<div align="center">
  <img src="05_media/02_prototype_alpha/side_view.jpg" alt="Prototype alpha top down view" style="max-width: 600px">
    <p><em>Figure 2: Prototype Alpha. Top-down view of the flat linear motor showing the slots & poles.</em></p>
</div>

The main conclusion from Prototype Alpha is that flat linear motors likely require laminated silicon steel armatures to produce force efficiently. In response, Prototype Beta shifts to an ironless tubular topology with the goal of quantifying force output and thermal performance.

> [!NOTE]  
> Prototype Alpha documentation and in-depth analysis can be found [here](02_motors/01_prototype_alpha/readme.md).

## Prototype Beta

<!-- Need to add a goal area for this like force target, input, asymptotic temperature, etc -->
