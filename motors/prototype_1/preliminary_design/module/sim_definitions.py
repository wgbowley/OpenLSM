"""
Filename: sim_definitions.py

Descriptions:
    Defines dataclasses that are used throughout
    the simulator and optimizer. This includes
    StaticEvaluation, DynamicSeries, DynamicEvaluation
    and MotorState.
"""


from __future__ import annotations

import csv

from pathlib import Path
from abc import ABC, abstractmethod
from dataclasses import dataclass, fields, field, InitVar

from pyfea import UnitError, Quantity as Q
from pyfea import (
    newton as N, ampere as A, second as S, henry as H, 
    weber as TM2, volt as V, kilogram as kg, meter as m, 
    ohm, kelvin as K, watt, dimensionless
)

 
_CREATED_THROUGH_CLS = object()


class SimData(ABC):
    """ Boundary data structures between modules """
    def validate_units(self) -> None:
        """ Generic validator that uses field metadata """
        for f in fields(self):
            required_unit = f.metadata.get('unit')
            if required_unit is None:
                continue
                
            attribute = getattr(self, f.name)
            
            if not isinstance(attribute, Q):
                msg = f"{f.name!r} must be a quantity, not {type(attribute)}"
                raise TypeError(msg)
            
            if attribute.unit != required_unit:
                msg = f"{f.name!r} must be {required_unit} not {attribute.unit}"
                raise UnitError(msg)

    @property 
    @abstractmethod
    def _name(self) -> str:
        """ Constructs a name based on attributes """

    def __post_init__(self) -> None: self.validate_units()
    def __repr__(self) -> str: return self._name


class SimStates(ABC):
    """Internal simulation state — construction restricted to .create()"""
    def __post_init__(self, _token: object) -> None:
        """ Fails direct user initialization """
        if _token is _CREATED_THROUGH_CLS:
            return

        name = self.__class__.__name__
        msg = f"{name} cannot be initialized directly. Use {name}.create()"
        raise PermissionError(msg)
    
    @property 
    @abstractmethod
    def _name(self) -> str:
        """ Constructs a name based on attributes """
        
    def __repr__(self) -> str: return self._name


@dataclass(frozen=True, slots=True, repr=False)
class StaticEvaluation(SimData):
    force_constant: Q           = field(metadata={'unit': N/A})
    resistance_atm_temp: Q      = field(metadata={'unit': ohm})
    secant_phase_inductance: Q  = field(metadata={'unit': H})
    magnet_flux: Q              = field(metadata={'unit': TM2})
    armature_mass: Q            = field(metadata={'unit': kg})
    slot_volume: Q              = field(metadata={'unit': m**3})
    armature_cost: Q            = field(metadata={'unit': dimensionless})
    segment_cost: Q             = field(metadata={'unit': dimensionless})
        
    @property
    def _name(self) -> str:
        """ Returns name as attributes """
        return (
            f"<StaticEvaluation("
            f"F/I: {self.force_constant:.3f}, "
            f"R_phase: {self.resistance_atm_temp:.3f}, "
            f"si_phase: {self.secant_phase_inductance:.3f}, "
            f"Φ_magnet: {self.magnet_flux}, "
            f"m_armature: {self.armature_mass:.3f}, "
            f"V_slot: {self.slot_volume:.3f}, "
            f"Armature Cost: ${self.armature_cost:.3f}, "
            f"Segment Cost: ${self.segment_cost:.3f})>"
        )
        

@dataclass(frozen=True, slots=True, repr=False)
class PathSegment(SimData):
    target_position: Q      = field(metadata={'unit': m})
    max_velocity: Q         = field(metadata={'unit': m/S})
    max_acceleration: Q     = field(metadata={'unit': m/S**2})
    time_out: Q             = field(metadata={'unit': S})

    @property
    def _name(self) -> str:
        """ Returns the name as attributes """
        return (
            f"<PathSegment("
            f"target_pos: {self.target_position:.3f}, "
            f"max_velocity: {self.max_velocity:.3f}, "
            f"max_acceleration: {self.max_acceleration:.3f}, "
            f"time_out: {self.time_out:.3f})>"
        )
    

@dataclass(slots=True, repr=False)
class DynamicSeries(SimStates):
    """ Time series results from dynamic analysis """
    time: Q
    set_points: Q
    k_m: Q
    displacement: Q
    velocity: Q
    d_current: Q
    q_current: Q
    force: Q
    total_power: Q
    p_loss: Q
    slot_temperature: Q
    pole_temperature: Q
    _token: InitVar[object] = field(default=None)

    @property
    def _name(self) -> str:
        """ Returns name as time series length """
        return f"<DynamicSeries(Entries: {len(self.time)})>"

    @classmethod
    def create(cls) -> DynamicSeries:
        """ Factory method to create set with units """
        return cls(
            []*S, []*m, []*(N/watt**0.5), []*m, []*m/S, []*A, 
            []*A, []*N, []*watt, []*watt, []*K, []*K, _CREATED_THROUGH_CLS
        )

    def record_step(
        self, state: MotorState, motor_constant: Q
    ) -> None:
        """ Records the new data by extracting fields from a MotorState instance """
        self.k_m.append(motor_constant)

        self.time.append(state.time)
        self.set_points.append(state.target_x)
        self.displacement.append(state.displacement)
        self.velocity.append(state.velocity)
        
        self.d_current.append(state.dq_currents[0])
        self.q_current.append(state.dq_currents[1])

        self.force.append(state.force)
        self.total_power.append(state.power)
        self.p_loss.append(state.power_loss)
        self.slot_temperature.append(state.slot_temperature)
        self.pole_temperature.append(state.pole_temperature)
        
    def to_csv(self, path: str | Path) -> None:
        """ Saves the time series data to a CSV file. """
        # Define the fields to export (mapping CSV headers to attribute names)
        fields = {
            "time (s)":            self.time,
            "set_points (m)":      self.set_points,
            "k_m (kg^1/2.s^-1/2)": self.k_m,
            "displacement (m)":    self.displacement,
            "velocity (m.s^-1)":   self.velocity,
            "d_current (A)":       self.d_current,
            "q_current (A)":       self.q_current,
            "force (N)":           self.force,
            "total_power (W)":     self.total_power,
            "p_loss (W)":          self.p_loss,
            "slot_temp (K)":       self.slot_temperature,
            "pole_temp (K)":       self.pole_temperature,
        }

        # Handle path expansion and resolution
        target_path = Path(path).expanduser().resolve()
        target_path.parent.mkdir(parents=True, exist_ok=True)

        # Extract magnitude data
        data_rows = zip(*(attr.value for attr in fields.values()))

        # Writing with 'utf-8' is still good practice, but now 
        # it won't crash even on 'cp1252' because everything is ASCII.
        with open(target_path, mode='w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(fields.keys())
            writer.writerows(data_rows)


@dataclass(frozen=True, slots=True, repr=False)
class DynamicEvaluation(SimData):
    motor_constant: Q           = field(metadata={'unit': N/watt**0.5})
    asymptotic_slot_temp: Q     = field(metadata={'unit': K})
    asymptotic_pole_temp: Q     = field(metadata={'unit': K})
    time_constant: Q            = field(metadata={'unit': H/ohm})
    weight: Q                   = field(metadata={'unit': kg})
    armature_cost: Q            = field(metadata={'unit': dimensionless})
    segment_cost: Q             = field(metadata={'unit': dimensionless})

    @property
    def _name(self) -> str:
        """ Returns name as attributes """
        return (
            f"<DynamicEvaluation("
            f"k_m: {self.motor_constant:.3f}, "
            f"asymptotic_slot_temp: {self.asymptotic_slot_temp:.3f}, "
            f"asymptotic_slot_temp: {self.asymptotic_pole_temp:.3f}, "
            f"T_constant: {self.time_constant}, "
            f"m_armature: {self.weight:.3f}, "
            f"Armature Cost: ${self.armature_cost:.3f}, "
            f"Segment Cost: ${self.segment_cost:.3f})>"
        )
    

@dataclass(slots=True)#, repr=False)
class MotorState(SimStates):
    """ Holds the motors state for an instance in time """
    # Time & control
    time: Q
    target_x: Q
    target_v: Q
    
    # Mechanical State
    force: Q
    velocity: Q
    acceleration: Q
    position_x: Q
    displacement: Q
    displace_angle: Q
    
    # Circuits
    inductance: Q
    resistance: Q
    
    # Electrical State
    dq_voltages: Q
    dq_induced: Q
    dq_currents: Q
    dq_flux_linkage: Q

    # Power & Flux State
    power: Q
    power_loss: Q
    total_flux: Q
    
    # Thermal State
    slot_temperature: Q
    pole_temperature: Q
    
    # State token
    _token: InitVar[object] = field(default=None)

    @property
    def _name(self) -> str:
        """ Returns name with time and position """
        return f"<MotorState(t={self.time}, x={self.displacement})>"

    @classmethod
    def create(cls, atm_temp: Q) -> MotorState:
        """ Factory method to create set with units """
        return cls(
            0*S, 0*m, 0*m/S,                # Time & Control
            0*N, 0*m/S, 0*m/S**2,           # Mech: F, v, a
            0*m, 0*m, 0*dimensionless,      # Mech: x, theta
            [0, 0]*H, 0*ohm,                # Circ: Ind_dq, Res
            [0, 0]*V, [0, 0]*V, [0, 0]*A,   # Elec: V_dq, V_ind, I_dq
            [0, 0]*TM2,                     # Elec: Flux_linkage
            0*watt, 0*watt, [0, 0, 0]*TM2,  # Power & Total Flux (New)
            atm_temp, atm_temp,             # Thermal 
            _CREATED_THROUGH_CLS
        )