"""
Filename: pd_pi_controller.py

Description:
    Cascaded PD-PI controller which controls the
    motor position via changing the driving voltage.
    
    This specific PD-PI controller is a proportional derivative position
    controller into a proportional integral current controller 
    with final outputs being a voltage between [0, supply_voltage]
"""

from math import sqrt

from modules.initial_setup import InitialConditions

from pyfea import Quantity as q, meter as M, ampere as A, volt as V, second as S
from pyfea.models.tubular_linear_motor.main import TubularLinearMotor

class CascadeController:
    """ Cascaded PD-PI controller for a linear motor """
    def __init__(
        self, motor: TubularLinearMotor, initial: InitialConditions
    ) -> None:
        """ Initializes the class and calculates loop terms """
        self.current_limit = motor.config.circuit.current_limit
        self.voltage_limit = motor.config.circuit.supply_voltage
        self.load = initial.armature_mass + motor.config.motion.load
        self.time_step = 0 * S

        # Defines controller constants
        freq = initial.resistance_atm_temp / initial.secant_phase_inductance

        self.freq = motor.config.numerical.de_solver_circuit_step * freq
        self.bandwidth = freq

        # Defines position loop variables
        self.loop_freq = self.freq / motor.config.numerical.pi_pd_solver_step
        self.damping_ratio = sqrt(2) / 2
        
        self.pos_kp = (self.load * self.loop_freq ** 2) / initial.force_constant
        self.pos_kd = (2 * self.load * self.damping_ratio * self.loop_freq ) / initial.force_constant
        self.target_position = 0 * M
        
        # Defines the current loop
        self.cur_kp = initial.secant_phase_inductance * self.bandwidth
        self.cur_ki = initial.resistance_atm_temp * self.bandwidth
        
        self.cur_summation = 0 * A * S
        self.cur_target = 0
        
    def set_target_position(self, position: q) -> None:
        """ Sets the new target position for the controller """
        self.target_position = position
    
    def step(
        self, position: q, velocity: q, target_velocity: q, current: q, force: q,
        d_induced: q
    ) -> q:
        """ Calculates the voltage to drive the motor to target position """
        self.cur_target = self._position_pd(
            position, velocity, target_velocity, force
        )
        return self._current_pi(current, d_induced)

    def sync_loop_time_step(self, time_step: q) -> None:
        """ Syncs up the Quasi-transient time_step and controller step """
        self.time_step = time_step
    
    def reset(self) -> None:
        """ Reset controller by wiping stateful variables """
        self.cur_summation = 0
        self.target_position = 0 * M
        self.cur_target = 0 
    
    def _position_pd(
        self, position: q, velocity: q, target_velocity: q, constant: q
    ) -> q:
        """ Calculates the current for the PI current controller """
        # self.pos_kp = (self.load * self.loop_freq ** 2) / (0.5 * constant)
        # self.pos_kd = (2 * self.load * self.damping_ratio *self.loop_freq) / (0.5 * constant)
        
        error = self.target_position - position
        proportional = self.pos_kp * error
        
        # Takes derivative on delta between target and actual velocity
        derivative = self.pos_kd * (target_velocity - velocity)
        current = proportional + derivative
        if abs(current) > self.current_limit:
            # Ensures that the directionality of the current is maintained
            return self.current_limit if current > 0 * A else - self.current_limit

        return current

    def _current_pi(self, current: q, d_induced: q) -> tuple[q, q]:
        """ Calculates the voltage for the motor q-axis; cancels d-axis voltage"""
        error = self.cur_target - current
        proportional = self.cur_kp * error
    
        # Predicted voltage with current integrator state
        tentative = proportional + self.cur_summation * self.cur_ki
        if abs(tentative) <= self.voltage_limit:
            self.cur_summation += error * self.time_step
            return (-d_induced.value, tentative.value) * tentative.unit

        # Camps and prevent windup
        windup = self.voltage_limit if tentative > 0 * V else - self.voltage_limit
        return (-d_induced.value, windup.value) * windup.unit