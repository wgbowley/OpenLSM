"""
Filename: model.py

Description:
    Model of tubular linear synchronous motor 
    for usage in transient simulation. 
"""

from typing import Any
from abc import ABC, abstractmethod
from builtins import float as f

from picounits import Quantity, Unit, strip_quantity, DynamicLoader
from picounits import TIME, VOLTAGE, NULLSET, LENGTH, DENSITY, CONDUCTIVITY, CURRENT, TEMPERATURE

from src.states import ArmatureState, PhaseState
from src.physics import derived

_, _ = ArmatureState, PhaseState

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
    """ Tubular linear synchronous motor model """
    def __init__(self, parameters: DynamicLoader, materials: DynamicLoader) -> None:
        """ Initializes the linear motor class """
        self.parameters = parameters
        self.materials = materials

        # Extract common variables and builds derived parameters
        self._extract_validate()
        self._construction()

    @property
    def _name(self):
        return super()._name

    def _construction(self) -> None:
        """ Constructs the tubular linear motor model """
        turns = derived.compute_turns(
            self.slot_axial_length, self.slot_thickness, self.slot_wire_diameter,
            self.slot_fill_factor
        )

        mean_slot_radius = (
            self.pole_thickness + self.armature_clearance + self.armature_thickness +
            self.slot_thickness / 2
        )
        resistance = derived.compute_resistance(
            turns, mean_slot_radius, self.slot_wire_diameter, 1/self.slot_conductivity
        )

        _ = resistance

    def _extract_validate(self) -> None:
        """ Extracts qualities from attribute tree and validates units """
        # Extract qualities from the numeric attribute tree
        numerics = self.parameters.numerics
        self.time_step = self.numericalize(numerics.time_step, TIME)
        self.line_voltage = self.numericalize(numerics.line_voltage, VOLTAGE)
        self.temperature = self.numericalize(numerics.temperature, TEMPERATURE)

        self.msg_frequency = self.numericalize(numerics.msg_freq, 1/TIME)

        # Extract qualities from the model attribute tree
        model = self.parameters.model
        self.pole_pairs = self.numericalize(model.pole_pairs, NULLSET)
        self.number_slots = self.numericalize(model.number_slots, NULLSET)

        # Extract qualities from the armature attribute tree
        armature = self.parameters.armature
        self.armature_slot_pitch = self.numericalize(armature.core.axial_slot_pitch, LENGTH)
        self.armature_clearance = self.numericalize(armature.core.radial_clearance, LENGTH)
        self.armature_thickness = self.numericalize(armature.core.radial_thickness, LENGTH)

        self.slot_wire_diameter = self.numericalize(armature.slots.wire_diameter, LENGTH)
        self.slot_fill_factor = self.numericalize(armature.slots.fill_factor, NULLSET)
        self.slot_axial_length = self.numericalize(armature.slots.axial_length, LENGTH)
        self.slot_thickness = self.numericalize(armature.slots.radial_thickness, LENGTH)

        # Extract qualities from the stator attribute tree
        stator = self.parameters.stator
        self.tube_thickness = self.numericalize(stator.tube.radial_thickness, LENGTH)
        self.pole_thickness = self.numericalize(stator.poles.radial_thickness, LENGTH)
        self.pole_axial_length = self.numericalize(stator.poles.axial_length, LENGTH)

        # Extracts slot qualities from material attribute tree
        self.slot_material = self.materials.copper
        self.slot_density = self.numericalize(self.slot_material.physical.density, DENSITY)
        self.slot_conductivity = self.linear_interpolate(
            self.slot_material.electrical.temperature_conductivity,
            numerics.temperature, CONDUCTIVITY
        )

        # Extracts pole qualities from material attribute tree
        self.pole_material = self.materials.NdFeB
        self.pole_coercivity = self._coercivity_at_temp(self.pole_material, numerics.temperature)

    @classmethod
    def linear_interpolate(cls, table: Quantity, value: Quantity, excepted: Unit) -> f:
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
                return cls.numericalize(result, excepted)

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
            return cls.numericalize(hc, CURRENT / LENGTH)

        return 0.0
