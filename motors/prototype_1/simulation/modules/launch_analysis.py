"""
Filename: launch_analysis.py

Description:
    Performs a Quasi-transient Electro-Magneto-Thermal-Mechanical
    analysis of a tubular linear motor launch between two z coords.
    
    NOTE:
    This scripts uses axial-symmetric coordinate system (Z-R) (Axial-Radial)
    due to modelling. Configuration for this launch can be found in 
    configuration.uiv
"""

from __future__ import annotations

from dataclasses import dataclass
from tabulate import tabulate
from operator import attrgetter

from modules.physics import *
from modules.initial_setup import InitialConditions
from modules.pd_pi_controller import CascadeController
from modules.suvat_feeding import SUVATFeeder

from pyfea import (
    Quantity as q, second as S, ampere as A, meter as M, kelvin as K,
    newton as N, joule as J, volt as V, weber, watt, ohm as Ω, henry as H,
    dimensionless, millimeter as mm
)

from pyfea.models.tubular_linear_motor.main import TubularLinearMotor
from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver
from pyfea.solver.femm.domains.thermostatic.solver import FEMMThermostaticSolver
from pyfea.solver.solver_outputs import (
    SolverOutputs, CircuitOptions, MagneticOptions, ThermalOptions
)


@dataclass(slots=True)
class SummaryLaunchResults:
    """ Summary results reports maximums """
    time: q
    maximum_currents: q
    maximum_velocity: q
    maximum_slot_temperature: q
    maximum_pole_temperature: q
    
    @property
    def _name(self) -> str:
        """ Returns its name as a table """
        table_data = [
            ["time", f"{self.time:.3f}"],
            ["Max Currents", f"{self.maximum_currents:.3f}"],
            ["Max Velocity", f"{self.maximum_velocity:.3f}"],
            ["Max Slot Temp", f"{self.maximum_slot_temperature:.3f}"],
            ["Max Pole Temp", f"{self.maximum_pole_temperature:.3f}"]
        ]

        table_str = tabulate(
            table_data, headers=["parameters", "value"], tablefmt="fancy_grid"
        )
    
        table_width = len(table_str.split('\n')[0])
        header_text = " Initial Conditions "
        
        return f"{header_text:=^{table_width}}\n{table_str}"
        
    def __repr__(self): return self._name
    

@dataclass(slots=True)
class LaunchResults:
    """ launch results holds the time series data """
    time: q
    d_current: q
    q_current: q
    force: q
    velocity: q
    displacement: q
    p_loss: q
    k_m: q
    slot_temperature: q
    pole_temperature: q
    set_points: q

    @classmethod
    def create(cls) -> LaunchResults:
        """ Factory method with units """
        return LaunchResults(
            []*S, []*A, []*A, []*N, []*(M/S), 
            []*M, []*watt, []*(N/watt**0.5), 
            []*K, []*K, []*M
        )

    def record_step(
        self, t: q, id: q, iq: q, f: q, v: q, d: q, 
        pl: q, km: q, st: q, pt: q, sp: q
    ) -> None:
        self.time.append(t)
        self.d_current.append(id)
        self.q_current.append(iq)
        self.force.append(f)
        self.velocity.append(v)
        self.displacement.append(d)
        self.p_loss.append(pl)
        self.k_m.append(km)
        self.slot_temperature.append(st)
        self.pole_temperature.append(pt)
        self.set_points.append(sp)


class _Analysis:
    """ Parent class for 'launch' defining low level primitives """
    def __init__(
        self, 
        motor: TubularLinearMotor,
        magnetic: FEMMMagnetostaticSolver,
        thermal: FEMMThermostaticSolver,
        initial_conditions: InitialConditions
    ) -> None:
        """ Initializes the class and defines dependencies """
        self.motor = motor
        self.initial_conditions = initial_conditions
        self.controller = CascadeController(motor, initial_conditions)
        self.magnetic = magnetic
        self.thermal = thermal
        
        # Quasi-transient loop variables
        self.results = LaunchResults.create()
        
        self.system_mass = (
            self.motor.config.motion.load + initial_conditions.armature_mass
        )

        # Magnetic requested outputs
        m_requested = SolverOutputs()
        m_requested.add_magnetic(self.motor.SLOT_ID, MagneticOptions.FORCE_LORENTZ)
        for phase in self.motor.PHASES:
            m_requested.add_circuit(phase, CircuitOptions.FLUX_LINKAGE)
            m_requested.add_circuit(phase, CircuitOptions.POWER)
            m_requested.add_circuit(phase, CircuitOptions.RESISTANCE)

        self.magnetic_outputs = m_requested
        
        # Thermal requested outputs
        t_requested = SolverOutputs()
        t_requested.add_thermal(self.motor.SLOT_ID, ThermalOptions.AVERAGE_TEMPERATURE)
        t_requested.add_thermal(self.motor.POLE_ID, ThermalOptions.AVERAGE_TEMPERATURE)
        
        self.thermal_outputs = t_requested
    
    def _solve_magnetic(
        self, displacement: q, angle: q, currents: q, 
        slot_temperature: q, pole_temperature
    ) -> tuple[q, q, q, q, q]:
        """ Solves the magnetic frame using the initialized solver"""    
        # Updates temperature of materials within the domain
        self.magnetic.update_temperature(
            self.motor.armature_slots_material, slot_temperature
        )
        self.magnetic.update_temperature(
            self.motor.stator_poles_material, pole_temperature
        )
        
        # Moves the armature by the vector (displacement, angle)
        moving_elements = [self.motor.SLOT_ID, self.motor.CORE_ID]
        self.magnetic.move_elements(moving_elements, displacement, angle)
        
        # Changes currents and add circuit results 
        for index, phase in enumerate(self.motor.PHASES):
            phase.current = currents[index]
            self.magnetic.update_current(phase)

        # Solves the magnetic problem and extracts results
        results = self.magnetic.solve(self.magnetic_outputs)
        
        power = 0 * watt
        resistance = 0 * Ω
        linkage = [] * weber
        for phase in self.motor.PHASES:
            power += abs(attrgetter(f"{phase.name}.power")(results))
            resistance += abs(attrgetter(f"{phase.name}.resistance")(results))
            linkage.append(attrgetter(f"{phase.name}.flux_linkage")(results))
        
        # Extracts force and turns to vector
        resistance /= len(self.motor.PHASES)
        force = attrgetter(f"element_{self.motor.SLOT_ID.value}.force_lorentz")(results)

        angle = 90 * dimensionless      # force[1].atan2(force[0]).to_degrees()
        force = force[1]                # Ignores radial forces

        return force, angle, resistance, power, linkage

    def _solve_thermal(self, power_loss: q, time_step: q) -> tuple[q, q]:
        """ Solves the thermal frame using the initialized solver """
        # Updates the volumetric heating based on power loss
        volumetric_heating = power_loss / self.initial_conditions.slot_volume
        self.thermal.update_heat_source(
            self.motor.armature_slots_material, volumetric_heating
        )
        
        # Solves the problem and than extracts pole and slot temperatures
        thermal_results = self.thermal.solve(self.thermal_outputs, time_step)
        pole = attrgetter(
            f"element_{self.motor.POLE_ID.value}.average_temperature"
        )(thermal_results)
        slot = attrgetter(
            f"element_{self.motor.SLOT_ID.value}.average_temperature"
        )(thermal_results)

        return slot, pole
        

class Launch(_Analysis):
    """ Runs the Quasi-transient of the linear motor between points """
    def run(self) -> tuple[SummaryLaunchResults, LaunchResults]:
        """ Performs a PD-PI controlled point to point analysis """
        inductance = self.initial_conditions.secant_phase_inductance
        resistance = self.initial_conditions.resistance_atm_temp
        magnet_flux = self.initial_conditions.magnet_flux
        system_mass = self.system_mass

        # Set the target for the controller
        displacement_target = self.motor.config.motion.axial_displacement
        feeder = SUVATFeeder(self.motor)
        
        # Fractional time constant stepping
        fraction = self.motor.config.numerical.de_solver_circuit_step
        time_step = inductance / (resistance) * fraction
        
        fraction = self.motor.config.numerical.thermal_step
        thermal_step = time_step * fraction
        
        self.controller.sync_loop_time_step(time_step)
        
        # Loop variables
        loop_time = 0 * S
        d_q_voltages = [0, 0] * V
        d_q_currents = [self.motor.config.numerical.initial_current, 0] * A
        a_b_c_currents = [0, 0, 0] * A
        t_flux = [0, 0, 0] * weber

        slot_temperature = self.motor.config.thermal.atmospheric_temperature
        pole_temperature = self.motor.config.thermal.atmospheric_temperature
        
        velocity = 0 * (M/S**1)
        displacement = 0 * M
        displacement_angle = 90 * dimensionless
        force = 0 * N
        power = 0.0 * watt

        target_x = 0 * M
        power_loss = 0 * watt
        prev_d_q_flux = [0, 0] * weber
        feeder.plan_path(displacement, displacement_target)
        while (
            0.05 * mm < abs(displacement - displacement_target) or
            0.05 * mm/S < abs(velocity)
        ):
            # Breaks the loop if it goes for too many time steps
            max = self.motor.config.numerical.de_solver_maximum_steps
            if int(loop_time.value / time_step.value) > max.value:
                break

            print (
                f"Step {loop_time/time_step:.0f}, Time: {loop_time:.3f}, "
                f"Power: {power:.3f}, Net Force: {force:.3f}, "
                f"velocity: {velocity:.3f}, dis {displacement:.3f}, tar: {target_x}, "
                f"dq_i: {d_q_currents}, dq_v: {d_q_voltages}"
            )

            # Updates the target position, velocity and feeds forwards velocity (phase leading)
            target_x, target_v = feeder.get_setpoint(loop_time)
            self.controller.set_target_position(target_x)
            
            # Calculates the d_q flux for the second loop as every loop is (n-1 behind)
            elec_angle = electrical_angle(displacement, self.motor.pole_pitch)
            a_b_current = inverse_park_transform(d_q_currents, elec_angle)
            a_b_c_currents = inverse_clarke_transform(a_b_current)
            
            a_b_frame = clarke_transform(t_flux[0], t_flux[1], t_flux[2])
            d_q_frame_flux = park_transform(a_b_frame, elec_angle)
            # Solves mechanical DE's through explicit euler method
            delta = velocity * time_step
            displacement += delta
            
            k_m = 0.0 * (N/watt ** 0.5)
            if power_loss > 0.0 * watt:
                k_m = force / power_loss ** 0.5

            self.results.record_step(
                loop_time, 
                d_q_currents[0],
                d_q_currents[1],
                force,
                velocity,
                displacement,
                power_loss,
                k_m,
                slot_temperature,
                pole_temperature,
                target_x
            )
            
            # Solves acceleration, velocity for the next frame
            acceleration = force / system_mass
            velocity += acceleration * time_step

            if velocity >= 0 * (M/S):
                displacement_angle = 90 * dimensionless 
            else:
                displacement_angle = 270 * dimensionless     
                
            # Solves magnetostatic frame using self.magnetic (solver)
            result = self._solve_magnetic(
                delta, displacement_angle, a_b_c_currents,
                slot_temperature, pole_temperature
            )

            # Extract solutions
            force, _ = result[0], result[1]
            resistance, _, t_flux = result[2], result[3], result[4]
            
            # Calculates the d-q flux linkage for current frame
            a_b_frame = clarke_transform(t_flux[0], t_flux[1], t_flux[2])
            d_q_frame_flux = park_transform(a_b_frame, elec_angle)
            
            d_q_ind = [inductance.value, inductance.value] * H
            if d_q_currents[0] > 0 * A and d_q_currents[1] > 0 * A:
                d_q_ind = [d_q_frame_flux[0]/d_q_currents[0], d_q_frame_flux[1]/d_q_currents[1]] * H
            
            current_magnitude = (d_q_currents[0]**2 + d_q_currents[1]**2) ** 0.5
            max_current = self.motor.config.circuit.current_limit
            
            if current_magnitude > max_current:
                # Emergency: scale back or halt
                scale = max_current / current_magnitude
                d_q_currents = [d_q_currents[0] * scale, d_q_currents[1] * scale]

            # Updates the system voltage via the pd-pi controller
            constant = self.initial_conditions.force_constant
            voltage = self.controller.step(
                displacement, velocity, target_v, d_q_currents[1], constant
            )
            d_q_voltages = [0, voltage.value] * voltage.unit

            # induced = (pi / self.motor.pole_pitch) * velocity * magnet_flux
            delta_flux_d = d_q_frame_flux[0] - prev_d_q_flux[0]
            delta_flux_q = d_q_frame_flux[1] - prev_d_q_flux[1]

            # d_induced_fea = delta_flux_d / time_step
            q_induced_fea = delta_flux_q / time_step

            induced = q_induced_fea  
            d_q_currents = rk_2nd_order_currents(
                d_q_currents,
                d_q_voltages,
                resistance,
                d_q_ind,
                induced,
                time_step
            )
            power = abs(d_q_voltages[0]*d_q_currents[0] + d_q_voltages[1]*d_q_currents[1])

            # Updates slot and pole temperature
            if loop_time.value % (thermal_step * time_step).value:
                power_loss = (d_q_currents[0]**2 + d_q_currents[1] ** 2) * resistance
                slot_temperature, pole_temperature = self._solve_thermal(
                    power_loss, time_step
                )

            prev_d_q_flux = d_q_frame_flux
            loop_time += time_step

        return None, self.results
                
            
            