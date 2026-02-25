"""
Filename: optimize.py

Description:
    Prototype 1 linear Motor optimization script refer to parameters.uiv
    for initial motor parameters.
    
    NOTE:
    Dependencies: pyfea, picounits, numpy, pymoo
    
    If this script fails to run. It most likely means that
    the pyfea (v0.1.0) api has been deprecated.
"""

import pickle
import numpy as np
from pathlib import Path
from pymoo.core.problem import ElementwiseProblem

from dataclasses import dataclass, fields

from module.sim_definitions import PathSegment
from module.initial_setup import static_evaluation
from module.dynamic_analysis import PointToPoint
from model.tubular import TubularLinearMotor

from pyfea import Quantity, dimensionless as D, millimeter as mm, second, kelvin, millisecond as ms, meter as m
from pyfea.solver.femm.domains.magnetostatic.solver import FEMMMagnetostaticSolver
from pyfea.solver.femm.domains.thermostatic.solver import FEMMThermostaticSolver

from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.optimize import minimize
from pymoo.operators.sampling.rnd import FloatRandomSampling

from pymoo.core.callback import Callback
from multiprocessing import Pool
from pymoo.parallelization import StarmapParallelization

# Defines configuration file and solver output path
BASE_DIR = Path(__file__).parent.parent.parent  
path_lib = BASE_DIR / "prototype_1/optimization/parameters.uiv"
solver_folder = BASE_DIR / "prototype_1/optimization/outputs"
motor_point_to_point = BASE_DIR / "prototype_1/optimization/path_data"

# Lists for indexed outputs:
pole_slot_ratios = [[2, 6] * D, [4, 12] * D, [6, 18] * D, [8, 24] * D]
pole_grades = ["N30", "N33", "N35", "N38", "N40", "N42", "N45", "N48", "N50", "N52"]
bare_conductor_diameters = [0.2, 0.224, 0.25, 0.315, 0.4, 0.5, 0.56, 0.71, 0.8, 1.0, 1.25, 1.5] * mm

# Configuration
segment = PathSegment(10 * mm, 200 * mm/second, 10000 * mm/second**2, 0.10 * second)
slot_temp_max = 393.15 * kelvin
pole_temp_max = 343.15 * kelvin
time_constant_max = 5 * ms
armature_length_max = 100 * mm


@dataclass
class InputsBounds:
    """ Holds the input parameters bounding domain """
    pole_slot_ratios:   tuple[Quantity, Quantity]   # (Indexed)
    poles_grade:        tuple[Quantity, Quantity]   # (Indexed)
    Wire_diameter:      tuple[Quantity, Quantity]   # (Indexed)
    slot_axial_length:  tuple[Quantity, Quantity]   # (continuous)
    slot_outer_radius:  tuple[Quantity, Quantity]   # (continuous)
    slot_axial_spacing: tuple[Quantity, Quantity]   # (continuous)
    
    def _collapse(self) -> tuple[list, list]:
        """ Collapses parameter domain into a list of lower, upper """
        lower, upper = [], []

        for field in fields(self):
            attribute = getattr(self, field.name)
            
            # Extracts quality and takes raw value
            lower.append(attribute[0].value)
            upper.append(attribute[1].value)

        return np.array(lower), np.array(upper)

    @property
    def collapse_lower(self) -> list:
        """ collapses parameter domain into lower bounds """
        return self._collapse()[0]
    
    @property
    def collapse_upper(self) -> list:
        """ collapses parameter domain into upper bounds """
        return self._collapse()[1]


class MyCheckpoint(Callback):
    def __init__(self, filepath="checkpoint.pkl"):
        super().__init__()
        self.filepath = filepath

    def notify(self, algorithm):
        # Save the entire algorithm state to a file
        with open(self.filepath, "wb") as f:
            pickle.dump(algorithm, f)
        print(f"  > Checkpoint saved to {self.filepath}")


def deletes_files(sim_name: str) -> None:
    """ Deletes simulation raw files from folder """
    thermal_raw = Path(solver_folder) / f"{sim_name}.feh"
    thermal_ans = Path(solver_folder) / f"{sim_name}.anh"
    thermal_raw.unlink()
    thermal_ans.unlink()
    
    magnetic_raw = Path(solver_folder) / f"{sim_name}.fem"
    magnetic_ans = Path(solver_folder) / f"{sim_name}.ans"
    magnetic_raw.unlink()
    magnetic_ans.unlink()


class OptimizationProblem(ElementwiseProblem):
    """ PYMOO evaluation wrapper for optimization. Evaluates one solution at a time """
    def __init__(self, bounds: InputsBounds, **kwargs) -> None:
        """ Entries optimization parameters into pymoo via super injection """
        super().__init__(
            n_var = 6,
            n_obj = 4,
            xl=bounds.collapse_lower,
            xu=bounds.collapse_upper,
            **kwargs
        )
        
    def _evaluate(self, x, out, *args, **kwargs) -> None:
        """ inherits the standard evaluate method from PYMOO problem """
        # Unpacks x variable into variables and creates a simulation name
        index1, index2, index3, axial_length, outer_radius, axial_spacing = x
        sim_name = f"sim_{index1}_{index2}_{axial_length}_{outer_radius}_{axial_spacing}"

        # Unpack & map indexed variables
        pole_slot       = pole_slot_ratios[int(round(index1))]
        pole_grade      = pole_grades[int(round(index2))]
        wire_diameter   = bare_conductor_diameters[int(round(index3))]
        
        # Creates motor and updates parameters
        TubularMotor = TubularLinearMotor(path_lib) 
        TubularMotor.update_parameters(
            pole_slot, pole_grade, wire_diameter, 
            axial_length * m, outer_radius * m, axial_spacing * m
        ) 
        
        # Initializes FEA solvers (magnetostatic and thermostatic)
        Magnetic = FEMMMagnetostaticSolver(solver_folder)
        Thermal = FEMMThermostaticSolver(solver_folder)
        
        f_penalty, g_penalty = [1e6] * 4, [1e6] * 4
        
        # Performs initializes evaluation and dynamic simulation for the motor
        try:
            static_results = static_evaluation(TubularMotor, Thermal, Magnetic, sim_name)
            
            # Quick force requirement check
            voltage = TubularMotor.params.circuit.supply_voltage
            current = voltage / static_results.resistance_atm_temp
            
            # Current limit (check)
            max_current = TubularMotor.params.circuit.current_limit
            current = max_current if max_current < current else current
            
            force = static_results.force_constant * current
            mass = static_results.armature_mass + TubularMotor.params.motion.load
            required_force = mass * segment.max_acceleration
            
            if force < required_force:
                deletes_files(sim_name)
                out["F"] = f_penalty
                out["G"] = [100, 100, 100, 100]
                return
            
            analysis = PointToPoint(TubularMotor, Thermal, Magnetic, static_results)

            csv_file = Path(motor_point_to_point) / f"{sim_name}.csv"
            _, results = analysis.run(segment, csv_file, True)
            
            # Extracts results
            motor_constant          = results.motor_constant.value
            asymptotic_slot_temp    = results.asymptotic_slot_temp.value
            asymptotic_pole_temp    = results.asymptotic_pole_temp.value
            time_constant           = results.time_constant.value
            # armature_weight         = results.weight.value
            armature_cost           = results.armature_cost.value
            pole_segment_cost       = results.segment_cost.value
            armature_length         = TubularMotor.armature_length.value
            
            deletes_files(sim_name)
            
            print(
                f" [OK] Km: {results.motor_constant:.4f} | "
                f"SlotT: {results.asymptotic_slot_temp:.1f} | "
                f"Tau: {results.time_constant:.2f} | "
                f"Cost: ${pole_segment_cost + armature_cost:.2f}"
            )
            
            # Objectives
            out["F"] = [
                -motor_constant,
                asymptotic_slot_temp,
                time_constant,
                pole_segment_cost + armature_cost
            ]
            
            # Constraints
            out["G"] = [
                asymptotic_pole_temp - pole_temp_max.value,
                asymptotic_slot_temp - slot_temp_max.value,
                time_constant - time_constant_max.value,
                armature_length - armature_length_max.value
            ]
                
        except Exception as err:
            deletes_files(sim_name)
            print(f"Simulation {sim_name!r} failed: {err}")
            out["F"] = f_penalty
            out["G"] = g_penalty


# Defines boundary of domains
bounds = InputsBounds(
    (0 * D,         (len(pole_slot_ratios) - 1) * D),
    (0 * D,         (len(pole_grades) - 1) * D),
    (0 * D,         (len(bare_conductor_diameters) -1) * D),
    (0.2 * mm,      10 * mm),
    (7.2 * mm,      14 * mm),
    (0.2 * mm,      10 * mm)
)

# Setups reference directions (NOTE: Not 100% on this)
ref_dirs = get_reference_directions("das-dennis", 4, n_partitions=4)

# Initializes algorithm
algorithm = NSGA3(
    pop_size=len(ref_dirs) + (4 - len(ref_dirs) % 4), # Multiples of 4 for parallel 
    ref_dirs=ref_dirs,
    sampling=FloatRandomSampling(),
)

if __name__ == "__main__":  
    pool = Pool(8)
    runner = StarmapParallelization(pool.starmap)
    problem = OptimizationProblem(bounds, elementwise_runner=runner)

    res = minimize(
        problem,
        algorithm,
        termination=('n_gen', 100),
        seed=1,
        save_history=True,
        verbose=True
    )

    print(" === Optimization Finished === ")
    print(f"Number of valid solutions found: {len(res.F)}")

    # Show the first 10 optimal motor designs
    for i in range(min(10, len(res.F))):
        print(f"\nSolution {i+1}:")
        print(f"Inputs: {res.X[i]}")
        print(f"Objectives: {res.F[i]}")