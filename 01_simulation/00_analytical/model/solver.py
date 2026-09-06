"""
Filename: solver.py

Description:
    Magnetic solver class for calculating the force of a 
    tubular linear motor using virtual work methods.
"""

from builtins import float as f
from math import pi, floor, sqrt

from picounits import Q, DynamicLoader, strip_quantity as validate
from picounits import LENGTH, VOLTAGE, CONDUCTIVITY, COERCIVITY, NULLSET, FORCE

from model.physics import field_equations
from model.physics import field_oriented_control


class Solver:
    """
    Magnetic Solver for linear tubular motor problem.
    Computes electromagnetic force using magnetic energy and virtual work methods.
    """

    def __init__(self, parameters: DynamicLoader) -> None:
        """ Initializes the solver class """
        self._extract_validate(parameters)
        self.permeability = 4 * pi * 10 ** -7

        # Phase currents
        self.i_pha = 0.0
        self.i_phb = 0.0
        self.i_phc = 0.0

        # Computes derived values from parameters
        self._compute_derived_values()

    def update_currents(self, angle: f) -> None:
        """ Updates phase currents based on time """
        ab_ref = field_oriented_control.inverse_park_transform(0, self.phase_peak_current, angle)
        a, b, c = field_oriented_control.inverse_clarke_transform(ab_ref)

        # Updates currents
        self.i_pha = a
        self.i_phb = b
        self.i_phc = c

    def compute_force(self, z_pos: Q) -> Q:
        """ Computes force using finite difference of the energy distribution """
        # Validates the unit of z_pos and removes units
        z_pos = validate(z_pos, LENGTH)

        electrical_angle = field_oriented_control.electrical_angle(z_pos, self.dipole_axial_length)
        self.update_currents(electrical_angle)

        # Calculates the energy one step forward and one backward
        pos = self._compute_energy_state(z_pos + self.der_step_size, self.int_step_size)
        neg = self._compute_energy_state(z_pos - self.der_step_size, self.int_step_size)

        force = - (pos - neg) / (2 * self.der_step_size)
        return force * FORCE

    def _compute_energy_state(self, translate: f, dz: f, epsilon: float = 1e-8, window: int = 5) -> f:
        """ Computes the co-energy interaction over the relevant domain. """
        positive = self._integrate_sample(dz, translate, epsilon, window)
        negative = self._integrate_sample(-dz, translate, epsilon, window)

        # Integrate cached values
        sum_energy = sum(positive) + sum(negative)
        energy = self.effective_area * self.permeability * sum_energy * dz
        return energy

    def _integrate_sample(self, dz: f, translate: f, epsilon: f, window: int) -> list[f]:
        """ Integrates the co-energy interaction in one direction """
        interactions = []
        z = 0.0
        below_epsilon = 0

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

        return interactions

    def _armature_field(self, pos: f = 0.0, translation: f = 0.0) -> f:
        """ Solves for the armature field at a specific z-position """
        phases = [self.i_pha, self.i_phb, self.i_phc]
        half_length = - self.number_slots * self.axial_slot_pitch / 2
        offset = self.armature_offset + half_length

        field_strength = 0.0
        for index in range(0, self.number_slots):
            # Selects the slots phase and polarity
            phase = phases[index % 3]
            polarity = -1 if index % 2 == 0 else 1

            # Computes the slot position
            slot_pos = offset + translation + self.axial_slot_pitch * index

            # Takes the sum of the fields at that position
            field_strength += field_equations.pole_field_strength(
                pos,
                slot_pos,
                phase,
                polarity * self.slot_turns,
                self.slot_axial_length,
                self.slot_mean_radius
            )

        return field_strength

    def _stator_dipole(self, index: int, pos: float, pole_pos: float) -> float:
        """ Single dipole within the motors stator """
        return (-1) ** index * (
            field_equations.dipole_field_strength(
                pos,
                pole_pos,
                self.coercivity,
                self.dipole_axial_length
            )
        )

    def _stator_field(self, pos: f = 0.0) -> f:
        """ Solves for the stator field at a specific z-position. """
        poles = int(self.tube_length / self.dipole_axial_length)
        half_poles = poles // 2

        # Include both positive and negative
        field_strength = 0
        for index in range(-half_poles + 1, half_poles):
            pole_pos = self.dipole_axial_length * index
            field_strength += self._stator_dipole(index, pos, pole_pos)

        return field_strength

    def _compute_derived_values(self) -> None:
        """ Compute derived values based on parameters """
        tube_outer_radius = self.dipole_radial_thickness + self.tube_radial_thickness
        core_inner_radius = tube_outer_radius + self.radial_clearance
        slot_inner_radius = core_inner_radius + self.core_radial_thickness

        # Calculates the mean radius & effective area
        self.slot_mean_radius = slot_inner_radius + self.slot_radial_thickness / 2
        self.effective_area = pi * self.slot_radial_thickness ** 2

        # Calculates the number of turns and inductance
        self.slot_turns = self._compute_turns()
        self.slot_inductance = self._compute_inductance()

        # Calculates the resistance and hence peak line current
        resistivity = 1 / self.conductivity
        self.slot_resistance = self._compute_resistance(resistivity)

        # Calculates the number of slots per phase & resistance
        self.phase_slots = self.number_slots // 3
        self.phase_resistance = self.phase_slots * self.slot_resistance
        self.phase_inductance = self.phase_slots * self.slot_inductance
        self.phase_voltage = self.line_voltage / sqrt(3)

        # Calculates the rms and peak currents for WYE
        self.phase_rms_current = self.phase_voltage / self.phase_resistance
        self.phase_peak_current = sqrt(2) * self.phase_rms_current

        # Calculate line-to-line resistance & inductance
        self.l2l_resistance = 2 * self.phase_resistance
        self.l2l_inductance = 2 * self.phase_inductance

        # Calculates the copper losses
        self.average_losses = 3 * self.phase_rms_current**2 * self.phase_resistance

        # Armature Offsets (z_0)
        self.armature_offset = self.dipole_axial_length / 2

    def _compute_turns(self) -> f:
        """ Computes the number of turns while according for the insulation & stacking. """
        slot_section = self.slot_axial_length * self.slot_radial_thickness
        wire_section = pi * (self.wire_diameter / 2) ** 2

        effective_area = slot_section * self.fill_factor
        return floor(effective_area / wire_section)

    def _compute_inductance(self) -> f:
        """ Calculates the coils self-inductance independent of mutual inductance between coils. """
        term = self.slot_turns ** 2 * self.permeability * self.effective_area
        return term / self.slot_axial_length

    def _compute_resistance(self, resistivity: f) -> f:
        """ Computes the slots resistance by using mean radius & conductor cross section. """
        turn_length = 2 * pi * self.slot_turns * self.slot_mean_radius
        cross_section = pi * (self.wire_diameter / 2) ** 2

        return resistivity * turn_length / cross_section

    def _extract_validate(self, parameters: DynamicLoader) -> None:
        """ Extracts qualities from attribute tree and validates units """
        # Numerical
        self.int_step_size = validate(parameters.numerics.solver.integration_step, LENGTH)
        self.der_step_size = validate(parameters.numerics.solver.derivative_step, LENGTH)
        self.line_voltage = validate(parameters.numerics.line_voltage, VOLTAGE)

        # Armature
        self.number_slots = validate(parameters.armature.number_slots, NULLSET)
        self.axial_slot_pitch = validate(parameters.armature.slots.axial_pitch, LENGTH)
        self.radial_clearance = validate(parameters.armature.radial_clearance, LENGTH)
        self.core_radial_thickness = validate(parameters.armature.core.radial_wall_thickness, LENGTH)

        # Armature Slots
        self.fill_factor = validate(parameters.armature.slots.material.fill_factor, NULLSET)
        self.wire_diameter = validate(parameters.armature.slots.material.wire_diameter, LENGTH)
        self.slot_axial_length = validate(parameters.armature.slots.axial_length, LENGTH)
        self.conductivity = validate(parameters.armature.slots.material.conductivity, CONDUCTIVITY)
        self.slot_radial_thickness = validate(parameters.armature.slots.radial_thickness, LENGTH)

        # Stator Tube
        self.tube_length = validate(parameters.stator.length, LENGTH)
        self.tube_radial_thickness = validate(parameters.stator.tube.radial_wall_thickness, LENGTH)

        # Stator Dipole
        self.coercivity = validate(parameters.stator.dipole.material.coercivity, COERCIVITY)
        self.dipole_axial_length = validate(parameters.stator.dipole.axial_length, LENGTH)
        self.dipole_radial_thickness = validate(parameters.stator.dipole.radial_thickness, LENGTH)
