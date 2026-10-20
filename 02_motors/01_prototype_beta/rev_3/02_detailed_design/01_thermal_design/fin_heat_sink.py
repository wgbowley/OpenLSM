"""
Filename: fin_heat_sink.py

Description:
    Designs the fins based on `scipython` modified 
    bessel functions with fixed temperature at base
    and thermal boundary at fin tips.
    
    NOTE:
    This is the initial thermal simulation.
    Goto 01_simulation for final thermal simulation.
"""


from math import exp, pi
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from scipy.special import k0, k1, i0, i1
from picounits import Q, Parser, resolve_derived

# Loads unit system, material library & parameters
ROOT_DIR = Path(__file__).resolve().parents[0]

resolve_derived()
parameters = Parser.open(ROOT_DIR / "parameters.uiv")
environmental = parameters.environmental

# Calculates the temperature at the fins
delta_fin = (
    parameters.source.copper_delta +
    parameters.source.potting_delta +
    parameters.source.base_delta
)

# Calculates the thermal expansion coefficient for an ideal gas
t_expansion = 1.0 / (delta_fin + parameters.environmental.temperature)


def compute_setup() -> tuple[Q, Q, Q, Q]:
    """ Computes the heat transfer, efficiency, number of fins and total wattage of the setup """
    radial_thickness = parameters.fin.outer_radius - parameters.fin.base_radius

    # Calculates the rayleigh number based on outer radius
    num = environmental.gravity * t_expansion * delta_fin * parameters.fin.axial_spacing ** 4
    den = environmental.k_viscosity * environmental.t_diffusivity * radial_thickness

    rayleigh = num / den

    # Computes the Elenbaas Correlation for natural convection
    nu_term = (1 - exp(-35/rayleigh)) ** (3/4)
    nu_natural = rayleigh / 24 * nu_term

    # Computes the heat transfer coefficient for the with the environment
    heat_transfer = nu_natural * environmental.k_conductivity / parameters.fin.axial_spacing

    # Cooling efficiency via modified bessel function
    beta = (heat_transfer / parameters.fin.k_conductivity / parameters.fin.axial_thickness) ** 0.5
    u0 = beta * parameters.fin.base_radius
    u1 = beta * parameters.fin.outer_radius

    first_term = k1(u0.stripped) * i1(u1.stripped) - i1(u0.stripped) * k1(u1.stripped)
    second_term = k0(u0.stripped) * i1(u1.stripped) + i0(u0.stripped) * k1(u1.stripped)

    bessel_function = first_term / second_term

    efficiency_term = beta * (parameters.fin.outer_radius**2 - parameters.fin.base_radius**2)
    efficiency = 2 * parameters.fin.base_radius / efficiency_term  * bessel_function

    # calculates the heat dissipation from a single fin
    fin_area = 2 * pi * (parameters.fin.outer_radius**2 - parameters.fin.base_radius**2)
    fin_heat_dissipation = efficiency * heat_transfer * fin_area * delta_fin

    # Calculates convective and radiate heat dissipation
    fins = parameters.fin.axial_length // (parameters.fin.axial_spacing + parameters.fin.axial_thickness)

    total_area = fin_area * fins
    radiate_term = (environmental.temperature + delta_fin)**4 - environmental.temperature**4
    radiate_dissipation = parameters.fin.emissivity * parameters.fin.stefan_boltzmann * total_area * radiate_term

    # Calculate exposed barrel height between fin attachments
    bare_height = parameters.fin.axial_length - (fins * parameters.fin.axial_thickness)
    bare_barrel_area = 2 * pi * parameters.fin.base_radius * bare_height
    bare_convection_dissipation = heat_transfer * bare_barrel_area * delta_fin

    # Weighting down radiate gains due self-absorption.
    heat_dissipation = (fins * fin_heat_dissipation) + bare_convection_dissipation + 0.5 * radiate_dissipation

    return heat_transfer, efficiency, fins, heat_dissipation


# Samples over the solution space using grid search
radial_solutions = []
for i in range(0, 100):
    axial_solutions = []

    parameters.fin.outer_radius += parameters.numerics.step_size.radius
    original = parameters.fin.axial_spacing

    for j in range(0, 100):
        parameters.fin.axial_spacing += parameters.numerics.step_size.axial_spacing

        _, _, _, result = compute_setup()
        axial_solutions.append((parameters.fin.axial_spacing, result))

    parameters.fin.axial_spacing = original
    radial_solutions.append((parameters.fin.outer_radius, axial_solutions))

    # Sample console output
    print(f"[{i + 1:3d}/100] radius: {parameters.fin.outer_radius:.4g}")


# Build arrays for plotting
radii = np.array([r.stripped for r, _ in radial_solutions])
spacings = np.array([s.stripped for s, _ in radial_solutions[0][1]])
heat_map = np.array([[r.stripped for _, r in axial] for _, axial in radial_solutions])

# Plot heat graph
fig, ax = plt.subplots(figsize=(9, 6))
im = ax.imshow(
    heat_map,
    origin="lower",
    aspect="auto",
    extent=[spacings[0], spacings[-1], radii[0], radii[-1]],
    interpolation="bilinear"
)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Heat dissipation (W)")

ax.set_xlabel("Axial spacing")
ax.set_ylabel("Outer radius")
ax.set_title("Heat dissipation over solution space")

fig.tight_layout()
plt.show()