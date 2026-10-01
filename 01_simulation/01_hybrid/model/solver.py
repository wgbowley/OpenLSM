"""
Filename: solver.py

Description:
    Magnetic solver class for calculating the force of a 
    tubular linear motor using virtual work methods.
"""


from dataclasses import dataclass

import numpy as np

from ifemm import Parser as iParser
from picounits import DynamicLoader, strip_quantity as validate
from picounits import LENGTH, VOLTAGE, CONDUCTIVITY, NULLSET


@dataclass(slots=True)
class Slot:
    """ A non-unit informed slot definition """
    turns: int
    current: float
    p1: tuple[float, float]
    p2: tuple[float, float]


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
        self.tube_radial_wall_thickness = validate(parameters.stator.tube.radial_wall_thickness, LENGTH)

        # Stator Dipole
        self.dipole_axial_length = validate(parameters.stator.dipole.axial_length, LENGTH)
        self.dipole_radial_thickness = validate(parameters.stator.dipole.radial_thickness, LENGTH)
