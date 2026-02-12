"""
Filename: initial_setup.py

Description:
    Performs a initial magnetostatic simulation of the 
    motor to get secant phase inductance, resistance @ atm_temp, 
    isolating magnet flux and phase alignment. 
    
    Solve thermal @ atm_temp to allow for transient thermal modelling
"""

from math import sin, sqrt, pi

from dataclasses import dataclass
from operator import attrgetter
from tabulate import tabulate

from pyfea.domain.units import Quantity, ampere as A, weber as W, ohm as Ω, henry as H
from pyfea.models.tubular_linear_motor.main import TubularLinearMotor
from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver
from pyfea.solver.femm.domains.thermostatic.solver import FEMMThermostaticSolver
from pyfea.solver.solver_outputs import (
    SolverOutputs, CircuitOptions, MagneticOptions, ThermalOptions
)

@dataclass(slots=True)
class InitialConditions:
    """ Holds initial conditions for quasi-transient loop """
    force_constant: Quantity
    resistance_atm_temp: Quantity
    secant_phase_inductance: Quantity
    magnet_flux: tuple[Quantity]
    armature_mass: Quantity
    slot_volume: Quantity
    
    @property
    def _name(self) -> str:
        """ Returns its name as a table """
        table_data = [
            ["Force Constant", f"{self.force_constant:.3f}"],
            ["Resistance @ Atm_Temp", f"{self.resistance_atm_temp:.3f}"],
            ["Secant Inductance", f"{self.secant_phase_inductance:.3f}"],
            ["Magnet Flux", f"{self.magnet_flux:.3f}"],
            ["Armature Mass", f"{self.armature_mass:.3f}"],
            ["Slot Volume", f"{self.slot_volume:.3f}"]
        ]

        table_str = tabulate(
            table_data, headers=["parameters", "value"], tablefmt="fancy_grid"
        )
    
        table_width = len(table_str.split('\n')[0])
        header_text = " Initial Conditions "
        
        return f"{header_text:=^{table_width}}\n{table_str}"
        
    
    def __repr__(self): return self._name

def pre_simulation_setup(
    motor: TubularLinearMotor,
    magnetic_solver: FEMMMagnetostaticSolver
) -> InitialConditions:
    """ Builds the initial condition dataclass via FEA / Analytical models """
    # Builds the simulation domain and sets phase currents to zero
    domain = motor.build_domain(magnetic_solver)
    magnetic_solver.setup(domain)
    
    for phase in motor.PHASES:
        phase.current = 0 * A
        magnetic_solver.update_current(phase)
    
    # Defines the required output for this simulation
    magnet_flux = SolverOutputs()
    for phase in motor.PHASES:
        magnet_flux.add_circuit(phase, CircuitOptions.FLUX_LINKAGE)
        
    # Solves and get average flux linkage across all phases
    magnet_results = magnetic_solver.solve(magnet_flux)
    magnet_flux = 0 * (W / A)
    for phase in motor.PHASES:
        magnet_flux += abs(
            attrgetter(f"{phase.name}.flux_linkage")(magnet_results)
        )
    magnet_flux /= len(motor.PHASES)    

    # Updates phases to a non-zero current position
    initial_current = motor.config.numerical.characterisation_current
    for index, phase in enumerate(motor.PHASES):
        trig_ratio = sin((5 * pi) / 6 + (4 * pi) / 3 * index)
        phase.current = initial_current * trig_ratio
        magnetic_solver.update_current(phase)
    
    # Defines the required output for this simulation
    general = SolverOutputs()
    general.add_magnetic(motor.SLOT_ID, MagneticOptions.FORCE_LORENTZ)
    general.add_magnetic(motor.SLOT_ID, MagneticOptions.VOLUME)
    general.add_magnetic(motor.CORE_ID, MagneticOptions.VOLUME)

    for phase in motor.PHASES:
        general.add_circuit(phase, CircuitOptions.FLUX_LINKAGE)
        general.add_circuit(phase, CircuitOptions.RESISTANCE)
        general.add_circuit(phase, CircuitOptions.CURRENT)

    # Solves and extracts resistance and inductance
    results = magnetic_solver.solve(general)
    secant_inductance, resistance = 0 * H, 0 * Ω
    for phase in motor.PHASES:
        flux = attrgetter(f"{phase.name}.flux_linkage")(results)
        current = attrgetter(f"{phase.name}.current")(results)
        resistance += attrgetter(f"{phase.name}.resistance")(results)
        secant_inductance += (abs(flux) - magnet_flux) / current
    
    resistance /= len(motor.PHASES)
    secant_inductance /= len(motor.PHASES)
    
    # Calculates the lorentz force constant
    current_rms = initial_current / sqrt(2)
    force = attrgetter(f"element_{motor.SLOT_ID.value}.force_lorentz")(results)
    force_constant = (force.magnitude * force.unit) / current_rms
    
    # Calculates armature mass
    slot_volume = attrgetter(f"element_{motor.SLOT_ID.value}.volume")(results)
    core_volume = attrgetter(f"element_{motor.SLOT_ID.value}.volume")(results)
    
    slot_density = motor.armature_slots_material.values().physical.density
    core_density = motor.armature_core_material.values().physical.density
    
    mass = slot_volume * slot_density + core_volume * core_density
    return InitialConditions(
        force_constant, resistance, secant_inductance, magnet_flux, mass, slot_volume
    )


def initial_state(tubular: TubularLinearMotor, solver: FEMMThermostaticSolver) -> None:
    """ Solves the initial state of the thermal problem without a heat source """
    domain = tubular.build_domain(solver)
    solver.setup(domain)

    outputs = SolverOutputs()
    outputs.add_thermal(tubular.POLE_ID, ThermalOptions.AVERAGE_TEMPERATURE)
    solver.solve(outputs)
