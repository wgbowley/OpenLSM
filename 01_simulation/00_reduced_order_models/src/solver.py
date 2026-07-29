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
        """ Computes force using finite difference of the energy distribution """
        # Calculates the energy one step forward and one backward
        pos = self._compute_energy_state(z_pos + dz, dz)
        neg = self._compute_energy_state(z_pos - dz, dz)

        return - (pos - neg) / (2 * dz)

    def _compute_energy_state(
        self, translate: f, dz: f, epsilon: float = 1e-8, window: int = 5
    ) -> f:
        """ Computes the co-energy interaction over the relevant domain. """
        positive = self._integrate_sample(dz, translate, epsilon, window)
        negative = self._integrate_sample(-dz, translate, epsilon, window)

        # Integrate cached values
        sum_energy = sum(positive) + sum(negative)
        energy = self.model.raw_effective_area * self.permeability * sum_energy * dz
        return energy

    def _integrate_sample(self, dz: f, translate: f, epsilon: f, window: int) -> list[f]:
        """ Integrates the co-energy interaction in one direction """
        interactions = []
        z = 0.0
        below_epsilon = 0

        while True:
            h_armature = self._armature_field(z, translate)
            h_stator = self._stator_field(z)
            interaction = h_armature * h_stator

            interactions.append(interaction)

            if abs(interaction) <= epsilon:
                below_epsilon += 1
            else:
                below_epsilon = 0

            # Breaks loop if below epsilon for more than window iterations
            if below_epsilon >= window: break

            # Moves along the z-axis
            z += dz

        return interactions

    def _armature_field(self, z_pos: f = 0.0, translation: f = 0.0) -> f:
        """ Solves for the armature field at a specific z-position """
        phases = [self.armature.phase_a, self.armature.phase_b, self.armature.phase_c]
        offset = - self.model.raw_armature_length / 2

        field_strength = 0
        for index in range(0, self.model.raw_number_slots):
            # Selects the slots phase and polarity
            phase = phases[index % 3]
            polarity = +1 if index % 2 == 0 else -1

            # Computes the slot position
            slot_start = offset + translation + self.model.raw_slot_pitch * index

            # Takes the sum of the fields at that position
            field_strength += dynamic.compute_slot_z_field_strength(
                z_pos,
                slot_start,
                phase.current,
                polarity * phase.turns,
                self.model.raw_slot_length,
                self.model.raw_slot_mean_radius
            )

        return field_strength

    def _stator_field(self, z_pos: f = 0.0) -> f:
        """ Solves for the stator field at a specific z-position. """
        poles = int(self.model.raw_stator_length / self.model.raw_pole_pitch)
        offset = - self.model.raw_stator_length / 2

        field_strength = 0
        for index in range(0, poles):
            # Computes the pole position
            pole_start = offset + self.model.raw_pole_pitch * index

            field_strength += (-1) ** index * dynamic.compute_stator_z_field_strength(
                z_pos,
                pole_start,
                self.model.raw_pole_pitch,
                self.model.raw_stator_field
            )

        return field_strength
