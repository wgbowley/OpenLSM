"""
Filename: dynamic_analysis.py

Description:
    Performs a Quasi-transient Electro-Magneto-Thermal-Mechanical
    analysis of a tubular linear motor launch between two z coords.
    
    NOTE:
    This scripts uses axial-symmetric coordinate system (Z-R) (Axial-Radial)
    due to modelling. Configuration for this launch can be found in 
    configuration.uiv
"""

from operator import attrgetter

from module.dynamic_physics import *
from module.suvat_feeding import SUVATFeeder
from module.pd_pi_controller import CascadeController
from model.tubular import TubularLinearMotor
from module.sim_definitions import (
    StaticEvaluation, DynamicSeries, DynamicEvaluation, MotorState, PathSegment
)

from pyfea import (
    Quantity as Q, ohm, weber, second, watt, newton, meter, dimensionless, henry, ampere, millisecond
)

from pyfea.solver.femm.domains.thermostatic.solver import FEMMThermostaticSolver
from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver
from pyfea.solver.solver_outputs import (
    SolverOutputs, CircuitOptions, MagneticOptions, ThermalOptions
)


class Analysis:
    """ Parent class for 'launch' defining thermal and magnetic simulations """
    def __init__(
        self, 
        motor: TubularLinearMotor,
        thermal: FEMMThermostaticSolver,
        magnetic: FEMMMagnetostaticSolver,
        static: StaticEvaluation
    ) -> None:
        """ Initializes the class and defines the dependencies """
        self.motor = motor
        self.static = static
        
        # Defines class dependencies
        self.controller = CascadeController(motor, static)
        self.thermal = thermal
        self.magnetic = magnetic
        
        # Defines system mass and required outputs
        self.system_mass = self.motor.params.motion.load + static.armature_mass
        self.series = DynamicSeries.create()
        
        # Magnetic requested outputs
        self.MagneticOutputs = SolverOutputs()
        self.MagneticOutputs.add_magnetic(motor.SLOT_ID, MagneticOptions.FORCE_LORENTZ)
        for phase in motor.PHASES:
            self.MagneticOutputs.add_circuit(
                phase, (CircuitOptions.FLUX_LINKAGE, CircuitOptions.RESISTANCE)
            )

        # Thermal requested outputs
        self.ThermalOutputs = SolverOutputs()
        self.ThermalOutputs.add_thermal(
            motor.POLE_ID, ThermalOptions.AVERAGE_TEMPERATURE
        )
        self.ThermalOutputs.add_thermal(
            motor.SLOT_ID, ThermalOptions.AVERAGE_TEMPERATURE
        )
        
    def _solve_magnetostatic(
        self, state: MotorState, currents: Q, delta: Q
    ) -> tuple[Q, Q, Q]:
        """ Solves a magnetostatic frame with updated frame values from state"""
        # Updates material properties to reflect temperature dependence materials
        self.magnetic.update_temperature(
            self.motor.armature_slots_material, state.slot_temperature
        )
        self.magnetic.update_temperature(
            self.motor.stator_poles_material, state.pole_temperature
        )
        
        # Moves the armature by the vector (displacement, angle)
        elements = [
            self.motor.SLOT_ID, 
            self.motor.CORE_ID, 
            self.motor.HEAT_SINK_ID, 
            self.motor.THERMAL_PASE_ID
        ]
        self.magnetic.move_elements(elements, delta, state.displace_angle)
        
        # Changes currents and solves the magnetic problem
        for index, phase in enumerate(self.motor.PHASES):
            phase.current = currents[index]
            self.magnetic.update_current(phase)
            
        # Solves the magnetic problem and extracts results
        results = self.magnetic.solve(self.MagneticOutputs)
        
        resistance = 0 * ohm
        linkage = [] * weber
        for phase in self.motor.PHASES:
            resistance += abs(attrgetter(f"{phase.name}.resistance")(results))
            linkage.append(attrgetter(f"{phase.name}.flux_linkage")(results))
        
        # Averages out any numerical errors across the phases
        resistance /= len(self.motor.PHASES)
    
        # Extracts axial force component 
        force = attrgetter(f"element_{self.motor.SLOT_ID.value}.force_lorentz")(results)[1]
        
        return force, resistance, linkage
    
    def _solve_thermal(self, power_loss: Q, time_step: Q = 0 * second) -> tuple[Q, Q]:
        """ Solves the thermal frame in either transient mode or asymptotic mode """
        # Calculates heating per m^3 and updates heat source
        volumetric = power_loss / self.static.slot_volume
        self.thermal.update_heat_source(self.motor.armature_slots_material, volumetric)
        
        # Solves the problem and than extracts pole and slot temperatures
        if time_step.value > 0:
            # Solves as a stateful transient model
            thermal_results = self.thermal.solve(self.ThermalOutputs, time_step)
        else:
            # Solves as a asymptotic model
            thermal_results = self.thermal.solve(self.ThermalOutputs)
            
        pole = attrgetter(
            f"element_{self.motor.POLE_ID.value}.average_temperature"
        )(thermal_results)
        slot = attrgetter(
            f"element_{self.motor.SLOT_ID.value}.average_temperature"
        )(thermal_results)
        
        return slot, pole


class PointToPoint(Analysis):
    """ Runs a Quasi-transient simulation of a linear motor between two z-points """
    def run(
        self, segments: tuple[PathSegment] | PathSegment, file_path: str, verbose: bool = False
    ) -> tuple[DynamicSeries, DynamicEvaluation]:
        tau = self.static.secant_phase_inductance / self.static.resistance_atm_temp

        # Time steps based on time constant of the linear motor at 1A
        main_time_step = 1/10 * tau
        
        # Setups controller via syncing loops
        self.controller.sync_loop_time_step(main_time_step)
        
        # Loop variables
        state = MotorState.create(self.motor.params.thermal.atmospheric_temperature)
        last_segment_finish = 0 * second
        
        # Initialize previous d-q flux for induced voltage calculation
        prev_dq_flux = [0, 0] * weber

        # Ensures segments is a list before beginning
        if not isinstance(segments, (tuple, list)):
            segments = [segments]
            
        for segment in segments:
            # Initializes the motion planner with this segments values
            feeder = SUVATFeeder(segment.max_velocity, segment.max_acceleration)
            feeder.plan_path(state.position_x, segment.target_position)
            
            not_meet_tolerance = True
            count = 0
            while not_meet_tolerance:
                # Breaks loop if maximum steps reached or segment time out is reached
                max_steps = self.motor.params.numerical.de_solver_maximum_steps
                if (state.time - last_segment_finish) / main_time_step > max_steps:
                    break
                elif (state.time - last_segment_finish) > segment.time_out:
                    break
                
                # Calculates the new targets for this step
                state.target_x, state.target_v = feeder.get_setpoint(
                    (state.time - last_segment_finish)
                )
                self.controller.set_target_position(state.target_x)
                
                if verbose:
                    print(
                        f"Step {state.time/main_time_step:.0f}, Time: {state.time:.3f}, "
                        f"Power: {state.power:.3f}, Force: {state.force:.3f}, "
                        f"Vel: {state.velocity:.3f}, Acc: {state.acceleration:.3f}, "
                        f"Pos: {state.position_x:.3f}, Tar: {state.target_x:.3f}, "
                        f"dq_v: {state.dq_voltages}, dq_i: {state.dq_currents}"
                    )

                # Calculates acceleration and velocity via explicit euler method
                normal_force = self.motor.params.motion.gravity * self.system_mass
                state.force -= normal_force * self.motor.params.motion.coefficient_friction
                
                state.acceleration = state.force / self.system_mass
                state.velocity += state.acceleration * main_time_step
                                    
                delta = state.velocity * main_time_step
                state.displacement += delta
                state.position_x = state.displacement
                
                motor_constant = 0.0 * (newton / watt ** 0.5)
                if state.power_loss > 0.0 * watt:
                    motor_constant = abs(state.force) / state.power_loss ** 0.5
                    
                # Records time series data for evaluation
                self.series.record_step(state, motor_constant)
                
                # Reflects directionality via (negative: left, positive: right)
                state.displace_angle = 270 * dimensionless
                if state.velocity >= 0 * meter / second:
                    state.displace_angle = 90 * dimensionless
                
                # Calculates the phase currents for this step
                elec_angle = electrical_angle(state.displacement, self.motor.pole_pitch)
                a_b_frame = inverse_park_transform(state.dq_currents, elec_angle)
                phase_currents = inverse_clarke_transform(a_b_frame)
                
                # Solves magnetostatic frame using inputted solver
                result = self._solve_magnetostatic(state, phase_currents, delta)
                state.force, state.resistance = result[0], result[1]
                state.total_flux = result[2]
                
                # NOTE: Due to numerical issues: Secant Inductance is USED
                # Calculates inductance
                # delta_flux = state.total_flux - self.static.magnet_flux
                # a_b_frame = clarke_transform(*delta_flux)
                # dq_delta_flux = park_transform(a_b_frame, elec_angle) 
                
                fallback = self.static.secant_phase_inductance
                state.inductance = [fallback.value, fallback.value] * henry

                # NOTE: Due to numerical issues: Secant Inductance is USED
                # if (
                #     abs(state.dq_currents[0]) > 0.1*ampere and 
                #     abs(state.dq_currents[0]) > 0.1*ampere
                # ):
                #     state.inductance = dq_delta_flux/state.dq_currents

                # Calculates the new d and q voltages to get to target
                state.dq_voltages = self.controller.step(state)
                
                # Solve for new currents using RK2 with induced voltages
                # This also updates prev_dq_flux internally
                state.dq_currents, prev_dq_flux = rk_2nd_order_currents(
                    self.motor,
                    state, 
                    prev_dq_flux,
                    main_time_step
                )
        
                # Calculates total and resistive loss power (scaled dq space)
                dq_power = state.dq_currents * state.dq_voltages
                state.power = 1.5 * (dq_power[0] + dq_power[1])
                state.power_loss = 1.5 * state.dq_currents.magnitude ** 2 * state.resistance
                
                # Updates slot and pole temperature
                thermal_count = self.motor.params.numerical.thermal_count.value
                if count == thermal_count:
                    slot, pole = self._solve_thermal(state.power_loss, thermal_count * main_time_step)
                    state.slot_temperature, state.pole_temperature = slot, pole
                    count = 0
                
                # Ensures that the armature actually has stopped.
                pos_tolerance = self.motor.params.motion.position_tolerance
                vel_tolerance = self.motor.params.motion.velocity_tolerance
                if (
                    pos_tolerance > abs(segment.target_position - state.position_x) and
                    vel_tolerance > abs(state.target_v - state.velocity)
                ):
                    not_meet_tolerance = False
                
                state.time += main_time_step
                count += 1
        
        # Saves dynamic data to csv format
        self.series.to_csv(file_path)
        
        # Calculates the asymptotic slot and pole temperature @ duty cycle
        time_on = 0.0 * second
        total_energy = 0.0 * (watt * second)
        for index, current in enumerate(self.series.q_current):
            if abs(current) > 0.1 * ampere:
                time_on += main_time_step
                total_energy += self.series.p_loss[index] * main_time_step
        
        average_power_loss = total_energy / state.time
        slot, pole = self._solve_thermal(average_power_loss)
        state.slot_temperature, state.pole_temperature = slot, pole

        # Calculates the average k_m when current is applied
        k_m = 0 * (newton / watt ** 0.5)
        count = 0
        for index, current in enumerate(self.series.q_current):
            if abs(current) > 0.1 * ampere:
                k_m += abs(self.series.k_m[index])
                count += 1

        k_m /= count

        # Given secant inductance usage throughout simulation t = L_s / R
        l_s = self.static.secant_phase_inductance / self.static.resistance_atm_temp
        return self.series, DynamicEvaluation(
            k_m, slot, pole, l_s, self.static.armature_mass, 
            self.static.armature_cost, self.static.segment_cost
        )