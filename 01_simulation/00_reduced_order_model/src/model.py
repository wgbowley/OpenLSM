"""
Filename: model.py

Description:
    Model of tubular linear synchronous motor 
    for usage in transient simulation. 
"""

from typing import Any
from abc import ABC, abstractmethod
from builtins import float as f


from math import pi
from picounits import Quantity, Unit, strip_quantity, DynamicLoader
from picounits import PERMEABILITY, INDUCTANCE, RESISTANCE, MASS, NULLSET
from picounits import LENGTH, CURRENT, VOLTAGE

from src.states import ArmatureState, PhaseState
from src.physics import derived


class SimulationModel(ABC):
    """ Abstract base class for simulation models """ 
    @classmethod
    def numericalize(cls, parameter: Quantity, expected: Unit) -> Any:
        """ Validates and strips a quantity into raw state """
        return strip_quantity(parameter, expected)

    @property
    @abstractmethod
    def _name(self) -> str:
        """ Constructs a name based on implementation attributes """

    def __repr__(self) -> str: return self._name
    def __str__(self) -> str: return str(self._name)


class TubularMotor(SimulationModel):
    """ 3 Phase Tubular linear Synchronous Motor """
    def __init__(self, parameters: DynamicLoader, materials: DynamicLoader) -> None:
        """ Initializes the model class """
        self.parameters = parameters
        self.materials = materials

        # Important variables
        self.permeability = 4 * pi * 10 ** -7 * PERMEABILITY
        self.temperature = self.parameters.numerics.temperature

        # Computes derived parameters
        self._compute_derived_phase_parameters()
        self.stator_field = self._coercivity_at_temp(self.materials.NdFeB, self.temperature)

        # Constructs the motor model
        self._construct_model()

        # Validates & Extracts Numerics Quantities
        self.raw_line_voltage = self.numericalize(self.parameters.numerics.line_voltage, VOLTAGE)

    def _construct_model(self) -> None:
        """ Validates, extracts and construct the model """
        params = self.parameters
        self.raw_number_slots = int(3 * self.numericalize(self.slots_per_phase, NULLSET))
        self.raw_pole_pitch = self.numericalize(params.stator.poles.axial_length, LENGTH)
        self.raw_stator_field = self.numericalize(self.stator_field, CURRENT/LENGTH)
        self.raw_stator_length = self.numericalize(params.stator.tube.length, LENGTH)

        self.raw_slot_pitch = self.numericalize(params.armature.core.axial_slot_pitch, LENGTH)
        self.raw_slot_length = self.numericalize(params.armature.slots.axial_length, LENGTH)
        self.raw_slot_mean_radius = self.numericalize(self.mean_radius, LENGTH)

        self.raw_effective_area = self.numericalize(self.effective_area, LENGTH ** 2)

        # Extracts the phase values
        raw_inductance = self.numericalize(self.phase_inductance, INDUCTANCE)
        raw_resistance = self.numericalize(self.phase_resistance, RESISTANCE)
        raw_mass = self.numericalize(self.phase_mass, MASS)
        raw_turns = self.numericalize(self.turns, NULLSET)

        # Calculates the armature length
        self.raw_armature_length = self.raw_slot_pitch * self.raw_number_slots


        # Constructs the phases & Armature state
        phase_a = PhaseState(raw_turns, raw_inductance, raw_resistance, raw_mass)
        phase_b = PhaseState(raw_turns, raw_inductance, raw_resistance, raw_mass)
        phase_c = PhaseState(raw_turns, raw_inductance, raw_resistance, raw_mass)

        self.armature = ArmatureState(phase_a, phase_b, phase_c)

    def _compute_derived_phase_parameters(self) -> None:
        """ Computes derived parameters for the motor phases """
        # Calculates the mean radius of turns within the slot
        params = self.parameters
        self.mean_radius = (
            params.stator.poles.radial_thickness +
            params.stator.tube.radial_thickness +
            params.armature.core.radial_clearance +
            params.armature.core.radial_thickness +
            params.armature.slots.radial_thickness / 2
        )

        # Calculates the number of turns and inductance
        self.turns = derived.compute_turns(
            params.armature.slots.axial_length,
            params.armature.slots.radial_thickness,
            params.armature.slots.wire_diameter,
            params.armature.slots.fill_factor
        )

        inductance = derived.compute_inductance(
            self.turns,
            params.armature.slots.axial_length,
            self.mean_radius,
            self.permeability
        )

        # Calculates resistivity at simulation temperature & slot resistance
        table = self.materials.copper.electrical.temperature_conductivity
        conductivity = self._linear_interpolate(table, self.temperature)
        resistivity = 1/conductivity

        resistance = derived.compute_resistance(
            self.turns,
            self.mean_radius,
            params.armature.slots.wire_diameter,
            resistivity
        )

        # Calculates slot volume and slot mass
        slot_volume = derived.compute_slot_volume(
            self.turns,
            params.armature.slots.wire_diameter,
            self.mean_radius
        )
        slot_mass = slot_volume * self.materials.copper.physical.density

        # Final inductance, resistance and mass per phase
        self.slots_per_phase = self.parameters.armature.core.number_slots / 3

        self.phase_inductance = self.slots_per_phase * inductance
        self.phase_resistance = self.slots_per_phase * resistance
        self.phase_mass = self.slots_per_phase * slot_mass

        # Calculates the interaction cross-sectional area
        self.effective_area = pi * self.mean_radius**2

    @classmethod
    def _linear_interpolate(cls, table: Quantity, value: Quantity) -> f:
        """ Linear interpolates a quantity list from a specific linked value """
        if value <= table[0][0]: return table[0][1]
        if value >= table[-1][0]: return table[-1][1]

        # finds specific interval
        for index in range(len(table) - 1):
            x0, x1 = table[index][0], table[index + 1][0]
            y0, y1 = table[index][1], table[index + 1][1]

            if x0 <= value <= x1:
                slope = (y1 - y0) / (x1 - x0)
                result = y0 + slope * (value - x0)
                return result

        msg = f"Value {value!r} could not be interpolated within table range {table!r}"
        raise ValueError(msg)

    @classmethod
    def _coercivity_at_temp(cls, material: DynamicLoader, temperature: Quantity) -> f:
        """ Calculates coercivity at a given temperature """
        ref_co = material.magnetic.coercivity
        beta_co = material.magnetic.beta_coercivity

        ref_temp = material.thermal.reference_temperature
        curie_point = material.thermal.curie_temperature

        if temperature < curie_point:
            hc = ref_co * (1 + (beta_co / 100) * (temperature - ref_temp))
            return hc

        return 0.0 * ref_co.unit

    @property
    def _name(self):
        """ Tubular motor representation """
        return "<TubularLinearSynchronousMotor(MagneticSolver)>"
