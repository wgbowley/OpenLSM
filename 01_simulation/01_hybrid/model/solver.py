"""
Filename: solver.py

Description:
    Magnetic solver class for calculating the force of a 
    tubular linear motor using virtual work methods.
"""


from builtins import float as f

import numpy as np

from ifemm import Parser as iParser
from picounits import DynamicLoader, strip_quantity as validate
from picounits import LENGTH, VOLTAGE, CONDUCTIVITY, NULLSET

from model.physics.field_equations import Slot


class Solver:
    """ Computes electromagnetic force using magnetic energy and virtual work methods. """
    def __init__(self, parameters: DynamicLoader, data: iParser) -> None:
        """ Initializes the solver class """
        self._extract_validate(parameters)

        # Finite element solution & environment permeability
        self.data = data
        self.permeability = 4 * np.pi * 10 ** -7

        # Phase currents
        self.i_pha = 0.0
        self.i_phb = 0.0
        self.i_phc = 0.0

        # Computes derived values from parameters
        self._compute_derived_values()

    def _armature_field(self, translation: f = 0.0) -> f:
        """ Constructs the armature field using slots """
        phases = [self.i_pha, self.i_phb, self.i_phc]
        half_length = - self.number_slots * self.slot_axial_length / 2
        offset = self.armature_offset + half_length

        slots = []
        for index in range(0, self.number_slots):
            # Selects the slots phase and polarity
            phase = phases[index % 3]
            polarity = -1 if index % 2 == 0 else 1

            # Computes the slot_pos & two reference points
            slot_pos = offset + translation + self.slot_axial_length * index
            p1 = (slot_pos, self.slot_inner_radius)
            p2 = (slot_pos + self.slot_axial_length, self.slot_outer_radius)

            # Appends the frozen slot
            slots.append(Slot(polarity * self.slot_turns, phase, p1, p2))

        return slots

    def _compute_derived_values(self) -> None:
        """ Compute derived values based on parameters """
        tube_outer_radius = self.dipole_radial_thickness + self.tube_radial_thickness
        core_inner_radius = tube_outer_radius + self.radial_clearance
        slot_inner_radius = core_inner_radius + self.core_radial_thickness
        slot_outer_radius = slot_inner_radius + self.slot_radial_thickness

        # Computes number of turns and saves slot size information
        self.slot_turns = self._compute_turns()

        self.slot_inner_radius = slot_inner_radius
        self.slot_outer_radius = slot_outer_radius

        # Armature Offsets (z_0)
        self.armature_offset = self.dipole_axial_length / 2

    def _compute_turns(self) -> f:
        """ Computes the number of turns while according for the insulation & stacking. """
        slot_section = self.slot_axial_length * self.slot_radial_thickness
        wire_section = np.pi * (self.wire_diameter / 2) ** 2

        effective_area = slot_section * self.fill_factor
        return int(np.floor(effective_area / wire_section))

    def _extract_validate(self, parameters: DynamicLoader) -> None:
        """ Extracts qualities from attribute tree and validates units """
        # Numerical Configuration
        self.line_voltage = validate(parameters.numerics.line_voltage, VOLTAGE)

        # Numerical Control
        self.integration_step = validate(parameters.numerics.solver.integration_step, LENGTH)
        self.derivative_step = validate(parameters.numerics.solver.derivative_step, LENGTH)

        # Armature Core
        self.number_slots = validate(parameters.armature.number_slots, NULLSET)
        self.radial_clearance = validate(parameters.armature.radial_clearance, LENGTH)
        self.core_radial_thickness = validate(parameters.armature.core.radial_wall_thickness, LENGTH)

        # Armature Slot
        self.slot_axial_pitch = validate(parameters.armature.slots.axial_pitch, LENGTH)
        self.slot_axial_length = validate(parameters.armature.slots.axial_length, LENGTH)
        self.slot_radial_thickness = validate(parameters.armature.slots.radial_thickness, LENGTH)

        self.wire_diameter = validate(parameters.armature.slots.material.wire_diameter, LENGTH)
        self.fill_factor = validate(parameters.armature.slots.material.fill_factor, NULLSET)
        self.conductivity = validate(parameters.armature.slots.material.conductivity, CONDUCTIVITY)

        # Stator Tube
        self.tube_radial_thickness = validate(parameters.stator.tube.radial_wall_thickness, LENGTH)

        # Stator Dipole
        self.dipole_axial_length = validate(parameters.stator.dipole.axial_length, LENGTH)
        self.dipole_radial_thickness = validate(parameters.stator.dipole.radial_thickness, LENGTH)
