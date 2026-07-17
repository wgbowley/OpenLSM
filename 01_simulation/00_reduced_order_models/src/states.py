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
            "<Armature("
            f"pos: {self.position:.3f} m, "
            f"vel: {self.position:.3f} m/s, "
            f"pha: {self.phase_a:.3f}, "
            f"phb: {self.phase_b:.3f}, "
            f"phc: {self.phase_c:.3f})>"
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
            f"vol: {self.voltage:.3f} V, "
            f"v_i: {self.induced_voltage:.3f} V, "
            f"Cur: {self.current:.3f} A, "
            f"Tur: {self.turns:.3f}, "
            f"Ind: {self.inductance:.3f} H, "
            f"Res: {self.resistance:.3f} Ω, "
            f"Mas: {self.mass:.3f} kg)>"
        )
