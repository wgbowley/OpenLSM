"""
Filename: optimize.py

Description:
    Prototype 1 linear Motor optimization script refer to parameters.uiv
    for initial motor parameters.
    
    NOTE:
    Dependencies: pyfea, picounits, numpy, pymoo
    
    If this script fails to run. It most likely means that
    the pyfea (v0.1.0) api has been deprecated.
    
    NOTE WINDOWS:
    $env:PYTHONUTF8=1
    [Console]::OutputEncoding = [System.Text.Encoding]::UTF8
    python -u file_path 2>&1 | tee output.log                                    
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
from pymoo.core.population import Population

from pymoo.core.callback import Callback
from multiprocessing import Pool
from pymoo.parallelization import StarmapParallelization

# Defines configuration file and solver output path
BASE_DIR = Path(__file__).parent.parent.parent  
path_lib = BASE_DIR / "prototype_1/optimization/parameters.uiv"
solver_folder = BASE_DIR / "prototype_1/optimization/outputs"
motor_point_to_point = BASE_DIR / "prototype_1/optimization/path_data"
checkpoint_path = BASE_DIR / "prototype_1/optimization/checkpoint.pkl"

# Lists for indexed outputs:
pole_slot_ratios = [[1, 6] * D, [2, 12] * D, [3, 18] * D, [4, 24] * D]
pole_grades = ["N30", "N33", "N35", "N38", "N40", "N42", "N45", "N48", "N50", "N52"]
bare_conductor_diameters = [0.2, 0.224, 0.25, 0.315, 0.4, 0.5, 0.56, 0.71, 0.8, 1.0, 1.25, 1.5] * mm

# Configuration
segment = PathSegment(50 * mm, 200 * mm/second, 10000 * mm/second**2, 0.5 * second)
slot_temp_max = 393.15 * kelvin
pole_temp_max = 343.15 * kelvin
time_constant_max = 10 * ms
armature_length_max = 200 * mm


@dataclass
class InputsBounds:
    """ Holds the input parameters bounding domain """
    pole_slot_ratios:           tuple[Quantity, Quantity]   # (Indexed)
    poles_grade:                tuple[Quantity, Quantity]   # (Indexed)
    Wire_diameter:              tuple[Quantity, Quantity]   # (Indexed)
    slot_axial_length:          tuple[Quantity, Quantity]   # (continuous)
    slot_radial_thickness:      tuple[Quantity, Quantity]   # (continuous)
    slot_axial_spacing:         tuple[Quantity, Quantity]   # (continuous)
    pole_radial_thickness:      tuple[Quantity, Quantity]   # (continuous)
    pole_axial_length:          tuple[Quantity, Quantity]   # (continuous)
    
    def _collapse(self) -> tuple[list, list]: 
        """ Collapses parameter domain into a list of lower, upper """
        lower, upper = [], []

        for field in fields(self):
            attribute = getattr(self, field.name)
            lower.append(attribute[0].value)
            upper.append(attribute[1].value)

        return np.array(lower), np.array(upper)

    @property
    def collapse_lower(self) -> list:
        return self._collapse()[0]
    
    @property
    def collapse_upper(self) -> list:
        return self._collapse()[1]


class MyCheckpoint(Callback):
    """
    Saves only the serialisable parts of the algorithm state
    """
    def __init__(self, filepath=checkpoint_path):
        super().__init__()
        self.filepath = Path(filepath)

    def notify(self, algorithm):
        checkpoint_data = {
            "n_gen": algorithm.n_gen,
            "pop_X": algorithm.pop.get("X"),
            "pop_F": algorithm.pop.get("F"),
            "pop_G": algorithm.pop.get("G"),
        }
        tmp_path = self.filepath.with_suffix(".tmp")
        with open(tmp_path, "wb") as f:
            pickle.dump(checkpoint_data, f)
        tmp_path.replace(self.filepath)
        print(f"  > Checkpoint saved (gen {algorithm.n_gen})", flush=True)


def deletes_files(sim_name: str) -> None:
    """ Deletes simulation raw files from folder """
    return
    # for ext in (".feh", ".anh", ".fem", ".ans"):
    #     (Path(solver_folder) / f"{sim_name}{ext}").unlink(missing_ok=True)


class OptimizationProblem(ElementwiseProblem):
    """ PYMOO evaluation wrapper for optimization. Evaluates one solution at a time """
    def __init__(self, bounds: InputsBounds, **kwargs) -> None:
        super().__init__(
            n_var=8,
            n_obj=6,
            n_constr=4,
            xl=bounds.collapse_lower,
            xu=bounds.collapse_upper,
            **kwargs
        )
        self.ref_force_constant   = 5
        self.ref_motor_constant   = 5
        self.ref_cost             = 50
        self.ref_position         = segment.target_position.value
        self.ref_time_constant    = time_constant_max.value
        self.ref_slot_temp        = slot_temp_max.value
        self.ref_pole_temp        = pole_temp_max.value
        self.ref_arm_length       = armature_length_max.value
        
    def _evaluate(self, x, out, *args, **kwargs) -> None:
        index1, index2, index3, axial_length, slot_thick, axial_spacing, pole_thick, pole_length = x

        pole_slot     = pole_slot_ratios[int(np.clip(np.round(index1), 0, len(pole_slot_ratios)-1))]
        pole_grade    = pole_grades[int(np.clip(np.round(index2), 0, len(pole_grades)-1))]
        wire_diameter = bare_conductor_diameters[int(np.clip(np.round(index3), 0, len(bare_conductor_diameters)-1))]
        
        sim_name = (
            f"sim_{index1:.1}_{index2:.1}_{index3:.1}_{axial_length:.6}_{slot_thick:.6}"
            f"_{axial_spacing:.6}_{pole_thick:.6}_{pole_length:.6}"
        )

        TubularMotor = TubularLinearMotor(path_lib) 
        TubularMotor.update_parameters(
            pole_slot, pole_grade, wire_diameter, 
            axial_length * m, slot_thick * m, axial_spacing * m, 
            pole_thick * m, pole_length * m
        ) 
        
        Magnetic = FEMMMagnetostaticSolver(solver_folder)
        Thermal  = FEMMThermostaticSolver(solver_folder)
        
        f_penalty = [10.0] * 6
        g_penalty = [10.0] * 4
        
        try:
            static_results = static_evaluation(TubularMotor, Thermal, Magnetic, sim_name)
            print(f"  [ok - static]  {static_results}", flush=True)

            time_constant = static_results.secant_phase_inductance / static_results.resistance_atm_temp
            
            voltage     = TubularMotor.params.circuit.supply_voltage
            current     = voltage / static_results.resistance_atm_temp
            max_current = TubularMotor.params.circuit.current_limit
            current     = max_current if max_current < current else current
            
            force         = static_results.force_constant * current
            mass          = static_results.armature_mass + TubularMotor.params.motion.load
            required_force = mass * segment.max_acceleration
            
            k_m_appox = force / (current * voltage).sqrt()
            
            # normalization of values
            force_norm = static_results.force_constant.value / self.ref_force_constant
            km_norm    = k_m_appox.value / self.ref_motor_constant
            tau_norm   = (
                static_results.secant_phase_inductance.value /
                static_results.resistance_atm_temp.value
            ) / self.ref_time_constant

            cost_norm = (
                static_results.armature_cost.value +
                static_results.segment_cost.value
            ) / self.ref_cost

            
            f_penalty = [
                -force_norm,      # obj 1
                -km_norm,         # obj 2
                10.0,             # slot temp (bad)
                10.0,             # position error (bad)
                tau_norm,         # time constant (scaled)
                cost_norm         # cost (scaled)
            ]

            g_penalty = [
                10.0,
                10.0,
                10.0,
                (TubularMotor.armature_length.value - self.ref_arm_length) / self.ref_arm_length
            ]
            
            if (force < required_force) or (tau_norm > 1.0):
                deletes_files(sim_name)
                out["F"] = f_penalty
                out["G"] = g_penalty
                return
            
            analysis = PointToPoint(TubularMotor, Thermal, Magnetic, static_results)
            csv_file = Path(motor_point_to_point) / f"{sim_name}.csv"
            dynamic, results = analysis.run(segment, csv_file)
            
            motor_constant       = results.motor_constant.value
            asymptotic_slot_temp = results.asymptotic_slot_temp.value
            asymptotic_pole_temp = results.asymptotic_pole_temp.value
            time_constant        = results.time_constant.value
            displacement         = dynamic.displacement[-1]
            armature_cost        = results.armature_cost.value
            pole_segment_cost    = results.segment_cost.value
            armature_length      = TubularMotor.armature_length.value
            
            target_different = segment.target_position - displacement
            deletes_files(sim_name)
            print(f"  [ok - dynamic]  {results} ", flush=True)
            
            out["F"] = [
                -force_norm,
                -motor_constant / self.ref_motor_constant,
                asymptotic_slot_temp / self.ref_slot_temp,
                abs(target_different.value) / self.ref_position,
                time_constant / self.ref_time_constant,
                (pole_segment_cost + armature_cost) / self.ref_cost
            ]

            out["G"] = [
                (asymptotic_pole_temp - self.ref_pole_temp) / self.ref_pole_temp,
                (asymptotic_slot_temp - self.ref_slot_temp) / self.ref_slot_temp,
                (time_constant - self.ref_time_constant) / self.ref_time_constant,
                (armature_length - self.ref_arm_length) / self.ref_arm_length
            ]

        except Exception as err:
            deletes_files(sim_name)
            print(f"Simulation {sim_name!r} failed: {err}", flush=True)
            out["F"] = f_penalty
            out["G"] = g_penalty


# Defines boundary of domains
bounds = InputsBounds(
    (0 * D,         (len(pole_slot_ratios) - 1) * D),
    (0 * D,         (len(pole_grades) - 1) * D),
    (0 * D,         (len(bare_conductor_diameters) - 1) * D),
    (0.2 * mm,      20 * mm),
    (7.2 * mm,      40 * mm),
    (0.2 * mm,      20 * mm),
    (2 * mm,        20 * mm),
    (2 * mm,        40 * mm)
)

ref_dirs = get_reference_directions("das-dennis", 6, n_partitions=5) 
pop_size = len(ref_dirs) + 1

if __name__ == "__main__":  
    pool   = Pool(1)
    runner = StarmapParallelization(pool.starmap)
    problem = OptimizationProblem(bounds, elementwise_runner=runner)

    initial_pop = None

    if checkpoint_path.exists():
        try:
            print(f"Resuming from checkpoint: {checkpoint_path}")
            with open(checkpoint_path, "rb") as f:
                data = pickle.load(f)
            print(f"  > Last completed gen: {data['n_gen']}", flush=True)
            # Warm-start from the last saved population
            initial_pop = data["pop_X"]
        except Exception as e:
            print(f"  > Checkpoint load failed ({e}), starting fresh...", flush=True)
            checkpoint_path.unlink(missing_ok=True)

    if initial_pop is None:
        print("Starting fresh...")
        random_samples = FloatRandomSampling().do(problem, pop_size).get("X")
        initial_pop = random_samples

    algorithm = NSGA3(
        pop_size=pop_size,
        ref_dirs=ref_dirs,
        sampling=initial_pop,
    )

    res = minimize(
        problem,
        algorithm,
        termination=('n_gen', 20),
        seed=1,
        save_history=True,
        verbose=True,
        callback=MyCheckpoint(checkpoint_path)
    )

    print(" === Optimization Finished === ", flush=True)
    print(f"Number of valid solutions found: {len(res.F)}", flush=True)

    for i in range(min(10, len(res.F))):
        print(f"\nSolution {i+1}:", flush=True)
        print(f"Inputs: {res.X[i]}", flush=True)
        print(f"Objectives: {res.F[i]}", flush=True)