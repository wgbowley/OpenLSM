"""
Filename: suvat_feeding.py

Description:
    Suvat feeder is drip feeds the pd-pi controller target
    position for it to meet during simulation. This allows
    for motion planning instead of directly targeting the
    final position. Also decrease overshot and oscillations.
"""


from pyfea import UnitError, Quantity as q, second, meter
from pyfea.models.tubular_linear_motor.main import TubularLinearMotor


class SUVATFeeder:
    """ Feeds a controller target position at target times """
    def __init__(self, motor: TubularLinearMotor) -> None:
        """ Initializes the feeder and adds key attributes """     
        # dp/dt, dv/dt maximums for the feeder to follow
        self.velocity_max = motor.config.motion.velocity_max
        self.acceleration_max = motor.config.motion.acceleration_max
    
    def plan_path(self, start_pos: q, end_pos: q) -> None:
        """ Calculates accelerations, velocities and time """
        if start_pos != meter or end_pos != meter:
            msg = f"Both start_pos and end_pos must be defined in meters"
            raise UnitError(msg)
        
        self.start_pos, self.dis = start_pos, abs(end_pos - start_pos)
        self.direction = 1 if (end_pos > start_pos) else -1

        # Calculates the time to reach v_max via v=u+at; assumes u=0
        self.t_accel = self.velocity_max / self.acceleration_max
        self.d_accel = 0.5 * self.acceleration_max * self.t_accel**2
        
        # Triangle profile (never hits maximum velocity)
        if self.d_accel * 2 > self.dis:
            self.d_accel = self.dis / 2
            self.t_accel = (2 * self.d_accel / self.acceleration_max) ** 0.5
            self.v_peak = self.acceleration_max * self.t_accel
            
            # Time at constant velocity
            self.t_cruise = 0 * second
        else:
            # Trapezoid profile (Hits maximum velocity; requires cruise period)
            self.v_peak = self.velocity_max
            
            # Calculates distance cruising than time to complete it
            self.t_cruise = (self.dis - 2*self.d_accel) / self.velocity_max
        
        self.total_time = 2 * self.t_accel + self.t_cruise
        
    def get_setpoint(self, loop_time: q) -> tuple[q, q]:
        """ Returns a position, velocity tuple for the current loop time """
        # Acceleration phase: s=1/2*at^2
        if loop_time < self.t_accel:
            v = self.acceleration_max * loop_time
            s = 0.5 * self.acceleration_max * loop_time ** 2
            
        # Cruise (constant v, a=0) s=u+v*t
        elif loop_time < (self.t_accel + self.t_cruise):
            t_offset = loop_time - self.t_accel
            v = self.v_peak
            s = self.d_accel + v * t_offset
            
        # Deceleration -> s=u+(ut+1/2*at^2)
        elif loop_time < self.total_time:
            t_offset = loop_time - self.t_accel - self.t_cruise
            v = self.v_peak - self.acceleration_max * t_offset
            
            s_decel = 0.5 * self.acceleration_max * t_offset ** 2
            s = (self.dis - self.d_accel) + (self.v_peak * t_offset - s_decel)
        
        # Move finished
        else:
            v = 0 * meter / second
            s = self.dis

        return self.start_pos + s * self.direction, v * self.direction