"""
Filename: sim_definitions.py

Descriptions:
    Defines dataclasses that are used throughout
    the simulator and optimizer. This includes
    StaticEvaluation, DynamicSeries, DynamicEvaluation
    and MotorState.
"""


from __future__ import annotations

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
            f"V_slot: {self.slot_volume:.3f})>"
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
            []*A, []*N, []*watt, []* watt, []*K, []*K, _CREATED_THROUGH_CLS
        )

    def record_step(
        self, t: Q, sp: Q, km: Q, d: Q, v: Q, id: Q, 
        iq: Q, f: Q, p: Q, pl: Q, st: Q, pt: Q
    ) -> None:
        """ Records the new data after a time step has occurred """
        self.time.append(t)
        self.set_points.append(sp)
        self.k_m.append(km)
        self.displacement.append(d)
        self.velocity.append(v)
        self.d_current.append(id)
        self.q_current.append(iq)
        self.force.append(f)
        self.total_power.append(p)
        self.p_loss.append(pl)
        self.slot_temperature.append(st)
        self.pole_temperature.append(pt)


@dataclass(frozen=True, slots=True, repr=False)
class DynamicEvaluation(SimData):
    motor_constant: Q       = field(metadata={'unit': N/watt**0.5})
    asymptotic_temp: Q      = field(metadata={'unit': K})
    force_ripple: Q         = field(metadata={'unit': dimensionless})
    time_constant: Q        = field(metadata={'unit': H/ohm})
    weight: Q               = field(metadata={'unit': kg})
    material_cost: Q        = field(metadata={'unit': dimensionless})

    @property
    def _name(self) -> str:
        """ Returns name as attributes """
        return (
            f"<DynamicEvaluation("
            f"k_m: {self.motor_constant:.3f}, "
            f"asymptotic_temp: {self.asymptotic_temp:.3f}, "
            f"Ripple: {self.force_ripple:.3f}, "
            f"T_constant: {self.time_constant}, "
            f"m_armature: {self.weight:.3f}, "
            f"mat_cost: {self.material_cost:.3f})>"
        )
    

@dataclass(slots=True, repr=False)
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
    displacement_angle: Q
    
    # Circuits
    inductance: Q
    resistance: Q
    
    # Electrical State
    dq_voltages: Q
    dq_induced: Q
    dq_currents: Q
    dq_flux_linkage: Q
    
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
            0*S, 0*m, 0*m/S,             # Time & Control
            0*N, 0*m/S, 0*m/S**2,        # Mech: F, v, a
            0 *m, 0*m, 0*dimensionless,  # Mech: x, theta
            [0, 0]*H, 0 * ohm            # Circ: Ind_dq, Res
            [0, 0]*V, [0, 0]*V, [0, 0]*A,# Elec: V_dq, I_dq
            [0, 0]*TM2,                  # Elec: Flux
            atm_temp*K, atm_temp*K,      # Thermal 
            _CREATED_THROUGH_CLS
        )