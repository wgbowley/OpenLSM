"""
Filename: initial_setup.py

Description:
    Performs a initial magnetostatic simulation of the
    motor to get values to define the pd_pi controller
    and loop variables.
    
    Also performs a thermostatic simulation @ atm_temp
    to initialize state for thermodynamics.
"""

from math import sin, pi

from module.sim_definitions import StaticEvaluation
from model.tubular import TubularLinearMotor

from pyfea import Quantity as Q, ampere as A, weber as TM2, henry as H, ohm


from pyfea.solver.femm.domains.thermostatic.solver import FEMMThermostaticSolver
from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver

from pyfea.solver.solver_outputs import (
    SolverOutputs, CircuitOptions, MagneticOptions, ThermalOptions
)

def _simulate_magnet_flux(
    motor: TubularLinearMotor,
    magnetic_solver: FEMMMagnetostaticSolver
) -> Q:
    """ Calculates magnet flux across phases """
    # Designs the phase currents at zero amperes
    for phase in motor.PHASES:
        phase.current = 0 * A
        magnetic_solver.update_current(phase)
        
    # Defines the required output for this simulation
    outputs = SolverOutputs()
    for phase in motor.PHASES: 
        outputs.add_circuit(phase, CircuitOptions.flux_linkage)
    
    # Solves and get average flux linkage across all phases
    magnet_results = magnetic_solver.solve(outputs)
    magnet_flux = [] * (TM2 / A)
    for phase in motor.PHASES:
        magnet_flux.append(magnet_results[phase].flux_linkage)

    return magnet_flux


def _simulate_magnetic_variables(
    motor: TubularLinearMotor,
    magnetic_solver: FEMMMagnetostaticSolver,
    magnet_flux: Q
) -> tuple[Q, Q, Q]:
    """ Calculates the resistance, inductance and force constant """
    characterisation_current = motor.params.numerical.characterisation_current
    
    # Updates phases to a non-zero current position
    for index, phase in enumerate(motor.PHASES):
        trig_ratio = sin(5 * pi / 6 + 4 * pi / 3 * index)
        phase.current = characterisation_current * trig_ratio
        magnetic_solver.update_current(phase)

    # Defines the required output for this simulation
    outputs = SolverOutputs()
    outputs.add_magnetic(motor.SLOT_ID, MagneticOptions.force_lorentz)

    for phase in motor.PHASES:
        outputs.add_circuit(phase, CircuitOptions.flux_linkage)
        outputs.add_circuit(phase, CircuitOptions.resistance)
        outputs.add_circuit(phase, CircuitOptions.current)    

    # Solves and extracts resistance and inductance
    results = magnetic_solver.solve(outputs)
    secant_inductance, resistance = 0 * H, 0 * ohm
    for index, phase in enumerate(motor.PHASES):
        flux = results[phase].flux_linkage
        current = results[phase].current
        resistance += results[phase].resistance

        secant_inductance += (abs(flux) - magnet_flux[index]) / current
    
    resistance /= len(motor.PHASES)
    secant_inductance /= len(motor.PHASES)
    
    # Calculates the lorentz force constant
    force = results[motor.SLOT_ID].force_lorentz
    force_constant = force[1] / characterisation_current
    
    return resistance, secant_inductance, force_constant


def initial_state(
    motor: TubularLinearMotor, solver: FEMMThermostaticSolver
) -> tuple[Q, Q, Q]:
    """ Solves the initial state of the thermal problem without a heat source """
    # Defines the required output for this simulation
    outputs = SolverOutputs()
    outputs.add_thermal(motor.SLOT_ID, ThermalOptions.volume)
    outputs.add_thermal(motor.CORE_ID, ThermalOptions.volume)
    outputs.add_thermal(motor.POLE_ID, ThermalOptions.volume)
    outputs.add_thermal(motor.HEAT_SINK_ID, ThermalOptions.volume)
    
    # Solves and extract parameters
    results = solver.solve(outputs)
    slot_volume = results[motor.SLOT_ID].volume
    core_volume = results[motor.CORE_ID].volume
    pole_volume = results[motor.POLE_ID].volume
    heat_volume = results[motor.HEAT_SINK_ID].volume

    # Material density
    slot_density = motor.armature_slots_material.values().physical.density
    core_density = motor.armature_core_material.values().physical.density
    pole_density = motor.stator_poles_material.values().physical.density
    heat_density = motor.heat_sink_material.values().physical.density
    
    # Material cost
    slot_cost = motor.params.material_cost.copper
    core_cost = motor.params.material_cost.nylon
    pole_cost = motor.params.material_cost.n52
    heat_cost = motor.params.material_cost.aluminum

    print(f"<Static volumes=(core={core_volume:.3f}, slot={slot_volume:.3f}, sink={heat_volume:.3f})")
        
    slot_mass, core_mass = slot_volume * slot_density, core_volume * core_density
    heat_mass = heat_volume * heat_density
    armature_mass = slot_mass + core_mass + heat_mass

    segment = pole_density * pole_volume * pole_cost
    armature = slot_mass * slot_cost + core_mass * core_cost + heat_cost * heat_mass
    return armature_mass, slot_volume, heat_volume, armature, segment


def static_evaluation(
    motor: TubularLinearMotor,
    thermal_solver: FEMMThermostaticSolver,
    magnetic_solver: FEMMMagnetostaticSolver,
    filename: str = "simulation_1"
) -> StaticEvaluation:
    """ Constructs 'StaticEvaluation' via simulation """
    # Builds the motor within the magnetostatic domain
    domain = motor.construct_domain(magnetic_solver)
    magnetic_solver.setup(domain, filename)
    
    # Calculates magnetic/circuit variables via FEA model
    magnet_flux = _simulate_magnet_flux(motor, magnetic_solver)
    resistance, inductance, force = _simulate_magnetic_variables(
        motor, magnetic_solver, magnet_flux
    )
    # Builds the motor within the thermostatic domain
    domain = motor.construct_domain(thermal_solver)
    thermal_solver.setup(domain, filename)
    armature_mass, slot_volume, heat_volume, armature, segment = initial_state(motor, thermal_solver)

    return StaticEvaluation(
        abs(force), abs(resistance), abs(inductance), magnet_flux, armature_mass, slot_volume, 
        heat_volume, armature, segment
    )