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

from model.physics.field_equations import kernel_limit, standard_helix, derivative_helix, biot_sum_integrand


class Solver:
    """ Computes electromagnetic force using magnetic energy and virtual work methods. """
    def __init__(self, parameters: DynamicLoader, data: iParser) -> None:
        """ Initializes the solver class """
        self._extract_validate(parameters)

        # Finite element solution, environment permeability & biot-savart term
        self.data = data
        self.permeability = 4 * np.pi * 10 ** -7
        self.mu0_over_4pi = self.permeability / (4 * np.pi)

        # Phase currents
        self.i_pha = 0.0
        self.i_phb = 0.0
        self.i_phc = 0.0

        # Stator field solution
        self.stator_r, self.stator_z, self.stator_bx, self.stator_by = data.b_field(512)
        self.b_mag_stator = np.sqrt(self.stator_bx**2 + self.stator_by**2)

        # Removing unit scaling from spacial axises
        self.stator_r *= data.length_scale
        self.stator_z *= data.length_scale

        # Computes derived values from parameters & slot kernel space
        self._compute_derived_values()
        self._construct_kernel()
        self._compute_kernel()

    def _compute_kernel(self) -> None:
        """ Computes the kernel solution based on slot geometry """
        effective_diameter = self.wire_diameter * (1 / self.fill_factor)

        # Calculates the number of sheets and turns per sheet.
        sheets = int(np.floor(self.slot_radial_thickness / effective_diameter))
        sheet_turns = int(np.floor(self.slot_axial_length / effective_diameter))

        # Creates the linear sample space for the sheets.
        total_samples = int(self.samples_per_turn * sheet_turns)
        t = np.linspace(0.0, 2 * np.pi * sheet_turns, total_samples)
        dt = t[1] - t[0]

        wires = []
        for sheet in range(1, sheets+1):
            # Calculates the new inner radius
            r_k = self.slot_inner_radius + effective_diameter * sheet

            # Computes the helix and its derivative
            wire = standard_helix(t, r_k, self.slot_axial_length, sheet_turns)
            d_wire = derivative_helix(t, r_k, self.slot_axial_length, sheet_turns, dt)

            # Appends wire & d_wire to wires set
            wires.append((wire, d_wire))

        # Constructs the kernel from integrand
        kernel = np.zeros_like(self.kernel_evaluation)
        for wire, d_wire in wires:
            kernel += biot_sum_integrand(self.kernel_evaluation, wire, d_wire)

        self.kernel = kernel
        print("finished")

    def _construct_kernel(self) -> None:
        """ Constructs the kernel size based on limit radius. """
        # Calculates the limit for the kernel
        nominal_radius = self.slot_inner_radius + self.slot_radial_thickness / 2
        limit = kernel_limit(self.slot_dropoff, self.slot_axial_length, nominal_radius)

        # Builds a mesh for the kernel with the same density as the FEM solution
        rx = round(1 / (self.stator_r[1] - self.stator_r[0]))
        rz = round(1 / (self.stator_z[1] - self.stator_z[0]))

        lin_r = np.linspace(0, limit, int(round(rx * limit)))
        lin_z = np.linspace(0, limit, int(round(rz * limit)))
        lin_y = np.linspace(0, limit, int(round(rz * limit)))

        # Creates the evaluation space for the kernel
        R, Y, Z = np.meshgrid(lin_r, lin_y, lin_z)
        r_eval = np.stack([R ,Y, Z], axis=-1)

        self.kernel_r = R
        self.kernel_y = Y
        self.kernel_z = Z

        self.kernel_evaluation = r_eval

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
        effective_turns = int(np.floor(effective_area / wire_section))

        if effective_turns <= 0:
            msg = f"A motor cannot have slots with negative turns: {effective_turns}"
            raise ValueError(msg)

        return effective_turns

    def _extract_validate(self, parameters: DynamicLoader) -> None:
        """ Extracts qualities from attribute tree and validates units """
        # Numerical Configuration
        self.line_voltage = validate(parameters.numerics.line_voltage, VOLTAGE)

        # Numerical Control
        self.slot_dropoff = validate(parameters.numerics.solver.slot_drop_off, NULLSET)
        self.samples_per_turn = validate(parameters.numerics.solver.samples_per_turn, NULLSET)

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
