"""
Filename: solver.py

Description:
    Solver class designed for solving 
    `TubularMotor` using virtual work methods.
"""

from builtins import float as f
from math import ceil, pi

from src.model import TubularMotor
from src.physics import dynamic


class MagneticSolver:
    """ 
    Magnetic field solver for `TubularMotor`.
    Computes electromagnetic force using magnetic energy and virtual work methods.
    """
    def __init__(self, model: TubularMotor) -> None:
        """ Initializes the solver class """
        self.model = model
        self.armature = self.model.armature

        # Important variables
        self.armature_length = self.model.raw_slot_pitch * self.model.raw_number_slots
        self.permeability = 4 * pi * 10 ** -7

    def compute_force(self, z_pos: f, dz: f) -> f:
        """  Computes force using finite difference of the energy distribution """
        # Calculates the energy one step forward and one backward
        pos = self._compute_energy_state(z_pos + dz, dz)
        neg = self._compute_energy_state(z_pos - dz, dz)

        return - (pos - neg) / (2 * dz)

    def _compute_energy_state(self, translate_pos: f, dz: f) -> f:
        """ Computes the energy distributed throughout the armature volume """
        energy = 0.0

        z_pos = - 1/2 * self.model.raw_stator_length
        for _ in range(ceil(self.model.raw_stator_length / dz)):
            h_armature = self._armature_field(z_pos - translate_pos)
            h_stator = self._stator_field(z_pos)

            # Calculates the energy using permeability
            energy += self.permeability * (h_armature + h_stator) ** 2 * dz
            z_pos += dz

        return (self.model.raw_effective_area/2) * energy

    def _armature_field(self, z_pos: f = 0.0) -> f:
        """ Solves for the armature field at a specific z-position """
        phases = [self.armature.phase_a, self.armature.phase_b, self.armature.phase_c]
        offset = - self.model.raw_slot_pitch * self.model.raw_number_slots / 2

        field_strength = 0
        for index in range(0, self.model.raw_number_slots):
            # Selects the slots phase and polarity
            phase = phases[index % 3]
            polarity = +1 if index % 2 == 0 else -1

            # Computes the slot position
            slot_start = offset + self.model.raw_slot_pitch * index
            slot_end = slot_start + self.model.raw_slot_length

            # Takes the sum of the fields at that position
            field_strength += dynamic.compute_slot_z_field_strength(
                z_pos,
                slot_end,
                phase.current,
                polarity * phase.turns,
                self.model.raw_slot_length,
                self.model.raw_slot_mean_radius
            )

        return field_strength

    def _stator_field(self, z_pos: f = 0.0, z_offset: f = 0.0) -> f:
        """ Solves for the stator field at a specific z-position """
        if abs(z_pos) > self.model.raw_stator_length / 2:
            return 0.0

        return dynamic.compute_stator_z_field_strength(
            z_pos,
            z_offset,
            self.model.raw_pole_pitch,
            self.model.raw_stator_field
        )
