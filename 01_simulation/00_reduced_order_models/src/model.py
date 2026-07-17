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

from src.states import ArmatureState, PhaseState
from src.physics import derived

_, _, _ = ArmatureState, PhaseState, derived


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

        # Constructs the model based on parameters & materials
        self._construction()

    @property
    def _name(self):
        return "LINEAR MOTOR"

    def _construction(self) -> None:
        """ Constructs the tubular linear motor model """

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
            return hc

        return 0.0 * ref_co.unit
