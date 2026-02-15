"""
Filename: pd_pi_controller.py

Description:
    Cascaded controller which controls the  motor position 
    via changing the driving voltage.
    
    This specific controller is a proportional derivative position
    controller into a proportional integral current controller 
    with final outputs being a voltage between [0, supply_voltage]
"""

from math import sqrt
from module.sim_definitions import StaticEvaluation, MotorState

from pyfea import UnitError, Quantity as Q
from pyfea import meter as m, ampere as A, volt as V, second as S
from pyfea.models.tubular_linear_motor.main import TubularLinearMotor


class CascadeController:
    """ Cascaded PD-PI motor controller """
    def __init__(self, motor: TubularLinearMotor, values: StaticEvaluation) -> None:
        """ Initializes the class and calculates the control variables """
        self.current_limit = motor.config.circuit.current_limit
        self.voltage_limit = motor.config.circuit.supply_voltage
        
        self.load = values.armature_mass + motor.config.motion.load
        self.time_step = 0 * S
        self.target_position = 0 * m
        
        # Defines controller frequency and bandwidth
        tau_frequency = values.resistance_atm_temp / values.secant_phase_inductance
        self.c_frequency = motor.config.numerical.de_solver_circuit_step * tau_frequency
        self.bandwidth = tau_frequency
        
        # Defines position loop variables
        self.l_frequency = self.c_frequency / motor.config.numerical.pi_pd_solver_step
        self.damping_ratio = sqrt(2) / 2
        
        f_constant = values.force_constant
        self.pos_kp = self.load * self.l_frequency ** 2 / f_constant
        self.pos_kd = 2 * self.load * self.damping_ratio * self.l_frequency / f_constant
        
        # Defines the current loop variables
        self.cur_kp = values.secant_phase_inductance * self.bandwidth
        self.cur_ki = values.resistance_atm_temp * self.bandwidth
        
        self.cur_integral = 0 * A*S
        self.cur_target = 0 * A

    def reset_variables(self) -> None:
        """ Reset controller by wiping stateful variables """
        self.target_position = 0 * m
        self.cur_integral, self.cur_target = 0 * A, 0 * A*S
            
    def set_target_position(self, position: Q) -> None:
        """ Sets the target position for the controller """
        if position.unit == m: 
            self.target_position = position
            return
            
        msg = f"Target position must be defined in meters, not {position.unit}"
        raise UnitError(msg)
    
    def sync_loop_time_step(self, time_step: Q) -> None:
        """ Syncs up the Quasi-transient time_step and controller step """
        if time_step == S: 
            self.time_step = time_step
            return
        
        msg = f"Time step must be defined in seconds, not {time_step.unit}"
        raise UnitError(msg)

    def step(self, state: MotorState) -> Q:
        """ Calculates the voltage to drive the motor to target position """
        self.cur_target = self._position_pd(state)
        return self._current_pi(self.cur_target)
        
    def _position_pd(self, state: MotorState) -> Q:
        """ Calculates the current for the PI controller from position/velocity """
        # Calculates the proportional of the position
        error = self.target_position - state.position_x
        proportional = self.pos_kp * error
        
        # Calculates the derivate of the velocity delta
        derivative = self.pos_kd * (state.target_v - state.velocity)
        current = proportional + derivative
        
        if abs(current) > self.current_limit:
            # Ensures that the directionality is maintained
            return self.current_limit if current > 0 * A else - self.current_limit
        
        return current
    
    def _current_pi(self, current: Q) -> Q:
        """ Calculates the voltage for motor q-axis from PD controller """
        # Calculates the proportional of the current
        error = self.cur_target - current
        proportional = self.cur_kp * error
        
        # Calculates the current integrator state
        tentative = proportional + self.cur_integral * self.cur_ki
        if abs(tentative) <= self.voltage_limit:
            self.cur_integral += error * self.time_step
            return tentative
        
        # Camps and prevent windup
        return self.voltage_limit if tentative > 0 * V else - self.voltage_limit