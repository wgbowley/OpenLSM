"""
Filename:

Description:
    Thermal solver class for calculating the total
    heat dissipation from an annulus heat sink.
"""

from math import floor, pi, log, exp
from builtins import float as f

from scipy.special import k0, k1, i0, i1

from picounits import Q, DynamicLoader, strip_quantity as validate
from picounits import TEMPERATURE, THERMAL_CONDUCTIVITY, LENGTH, TIME, NULLSET, POWER

class Solver:
    """ Computes the total heat dissipation for a annulus heat sink """
    def __init__(self, parameters: DynamicLoader) -> None:
        """ Initializes the solver class """
        self._extract_validate(parameters)
        self.stefan_boltzmann = 5.67 * 10 ** - 8

        # Computes derived values from parameters
        self._compute_derived_values()

    def compute_performance(self) -> tuple[Q, Q, Q]:
        """ Computes the thermal performance map of a series of heat-sinks """
        x_axial, y_radial, z_dissipation = [], [], []

        thickness = self.radial_size
        for _ in range(self.samples):
            y_radial.append(thickness)
            thickness += self.radial_size

        axial = self.axial_size
        for i in range(self.samples):
            x_axial.append(axial)

            row = []
            thickness = self.radial_size
            for _ in range(self.samples):
                dissipation = self._try_permutation(axial, thickness)
                row.append(dissipation)

                thickness += self.radial_size

            z_dissipation.append(row)
            axial += self.axial_size

            # Sample console output
            print(f"[{i + 1:3d}/100] spacing: {axial * LENGTH}")

        # Returns (m), (m), (W=kg*m^2*s^-3)
        return x_axial, y_radial, list(zip(*z_dissipation))

    def compute_permutation(self, axial_spacing: Q, radial_thickness: Q) -> Q:
        """ Computes a specific geometric permutation based on axial space and radial thickness """
        axial = validate(axial_spacing, LENGTH)
        radial = validate(radial_thickness, LENGTH)

        # Computes the permutation and returns
        dissipation = self._try_permutation(axial, radial)
        return dissipation * POWER

    def _try_permutation(self, axial_spacing: f, radial_thickness: f) -> f:
        """ Computes the thermal performance of the axial space and thickness permutation """
        fin_outer_radius = self.sink_outer_radius + radial_thickness

        # Computes the rayleigh number based on spacing and thickness
        rayleigh = self.gravity * self.t_expansion * self.fin_delta * axial_spacing ** 4
        rayleigh /= self.viscosity * self.diffusivity * radial_thickness

        # Computes the elenbaas correlation for natural convection and heat transfer coefficient
        nu_natural = rayleigh / 24 * (1 - exp(-35 / rayleigh)) ** (3/4)
        heat_transfer = nu_natural * self.env_conduct / axial_spacing

        # Cooling efficiency via modified bessel function
        beta = (heat_transfer / self.sink_conductivity / self.fin_axial_thickness) ** (1/2)
        u0 = beta * self.sink_outer_radius
        u1 = beta * fin_outer_radius

        first_term = k1(u0) * i1(u1) - i1(u0) * k1(u1)
        second_term = k0(u0) * i1(u1) + i0(u0) * k1(u1)

        bessel_function = first_term / second_term

        efficiency = 1 / (beta * (fin_outer_radius ** 2 - self.sink_outer_radius ** 2))
        efficiency *= 2 * self.sink_outer_radius * bessel_function

        # Calculates the heat dissipation from a single fin
        fin_area = 4 * pi * (fin_outer_radius ** 2 - self.sink_outer_radius ** 2)
        per_fin_heat_dissipation = efficiency * heat_transfer * fin_area * self.fin_delta

        # Calculates the number of fins and total area
        fins = floor(self.sink_axial_length // (axial_spacing + self.fin_axial_thickness))

        total_fin_area = fin_area * fins
        bare_barrel_area = 2 * pi * self.sink_outer_radius * (self.sink_axial_length - fins * self.fin_axial_thickness)
        total_area = bare_barrel_area + total_fin_area

        # Calculates the radiate heat dissipation and convection dissipation
        radiate = (self.env_temp + self.fin_delta) ** 4 - self.env_temp ** 4
        radiate *= self.sink_emissivity * self.stefan_boltzmann * total_area

        convection = heat_transfer * bare_barrel_area * self.fin_delta + fins * per_fin_heat_dissipation

        # Total heat dissipation
        total_dissipation = convection + self.fin_self_absorption * radiate
        return total_dissipation

    def _compute_derived_values(self) -> None:
        """ Compute derived values based on parameters """
        potting_outer_radius = self.potting_inner_radius + self.potting_radial_thickness
        self.sink_outer_radius = potting_outer_radius + self.sink_radial_thickness

        # Calculates the temperature drop across the potting material
        r_pot = log(potting_outer_radius / self.potting_inner_radius)
        r_pot /= 2 * pi * self.potting_conductivity * self.sink_axial_length

        self.t_pot_drop = - r_pot * self.copper_loss

        # Calculates the temperature drop across the annulus heat-sink base
        r_sink = log(self.sink_outer_radius / potting_outer_radius)
        r_sink /= 2 * pi * self.sink_conductivity * self.sink_axial_length

        self.t_sink_drop = - r_sink * self.copper_loss
        self.fin_delta = self.copper_delta + self.t_pot_drop + self.t_sink_drop

        # Calculates the thermal expansion coefficient (Assumes ideal gas)
        self.t_expansion = 1.0 / (self.fin_delta + self.env_temp)

    def _extract_validate(self, parameters: DynamicLoader) -> None:
        """ Extracts qualities from attribute tree and validates units """
        # Numerics
        self.samples = validate(parameters.numerics.samples, NULLSET)
        self.axial_size = validate(parameters.numerics.step_size.axial_spacing, LENGTH)
        self.radial_size = validate(parameters.numerics.step_size.radial_thickness, LENGTH)

        self.copper_loss = validate(parameters.numerics.source.loss, POWER)
        self.copper_delta = validate(parameters.numerics.source.delta, TEMPERATURE)

        # Environmental variables
        self.gravity = validate(parameters.environmental.gravity, LENGTH * TIME ** -2)
        self.env_temp = validate(parameters.environmental.temperature, TEMPERATURE)
        self.viscosity = validate(parameters.environmental.viscosity, LENGTH ** 2 * TIME ** -1)
        self.diffusivity = validate(parameters.environmental.diffusivity, LENGTH ** 2 * TIME ** -1)
        self.env_conduct = validate(parameters.environmental.conductivity, THERMAL_CONDUCTIVITY)

        # Potting
        self.potting_inner_radius = validate(parameters.potting.inner_radius, LENGTH)
        self.potting_radial_thickness = validate(parameters.potting.radial_thickness, LENGTH)
        self.potting_conductivity = validate(parameters.potting.material.conductivity, THERMAL_CONDUCTIVITY)

        # Heat Sink
        self.sink_axial_length = validate(parameters.sink.axial_length, LENGTH)
        self.sink_radial_thickness = validate(parameters.sink.radial_thickness, LENGTH)
        self.fin_axial_thickness = validate(parameters.sink.fin.axial_thickness, LENGTH)
        self.fin_self_absorption = validate(parameters.sink.fin.self_absorption, NULLSET)
        self.sink_emissivity = validate(parameters.sink.material.emissivity, NULLSET)
        self.sink_conductivity = validate(parameters.sink.material.conductivity, THERMAL_CONDUCTIVITY)
