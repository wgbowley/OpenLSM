"""
Filename: states.py

Descriptions:
    Motor states dataclasses for usage
    in tracking transient values.
"""

from __future__ import annotations

from builtins import float as f

from abc import ABC, abstractmethod
from dataclasses import dataclass

from picounits import VOLTAGE, CURRENT, INDUCTANCE, RESISTANCE, NULLSET, MASS, LENGTH, TIME


class SimulationState(ABC):
    """ Abstract base class for simulation states """
    @property
    @abstractmethod
    def _name(self) -> str:
        """ Constructs a name based on implementation attributes """

    def __repr__(self) -> str: return self._name
    def __str__(self) -> str: return str(self._name)


@dataclass(slots=True)
class ArmatureState(SimulationState):
    """ Holds variables related to the armature """
    phase_a: PhaseState
    phase_b: PhaseState
    phase_c: PhaseState
    velocity: f = 0.0
    position: f = 0.0

    @property
    def _name(self):
        """ Returns name as attributes """
        return (
            "<Armature(\n"
            f"pos: {self.position * LENGTH:.3f}, \n"
            f"vel: {self.position * (LENGTH / TIME):.3f}, \n"
            f"pha: {self.phase_a}, \n"
            f"phb: {self.phase_b}, \n"
            f"phc: {self.phase_c}\n)>"
        )


@dataclass(slots=True)
class PhaseState(SimulationState):
    """ Holds variables related to a single phase """
    turns: f
    inductance: f
    resistance: f
    mass: f
    voltage: f = 0.0
    induced_voltage: f = 0.0
    current: f = 0.0

    @property
    def _name(self):
        """ Returns name as attributes """
        return (
            "<Phase("
            f"vol: {self.voltage * VOLTAGE:.3f}, "
            f"v_i: {self.induced_voltage * VOLTAGE:.3f}, "
            f"Cur: {self.current * CURRENT:.3f}, "
            f"Tur: {self.turns * NULLSET:.3f}, "
            f"Ind: {self.inductance * INDUCTANCE:.3f}, "
            f"Res: {self.resistance * RESISTANCE:.3f}, "
            f"Mas: {self.mass * MASS:.3f})>"
        )
