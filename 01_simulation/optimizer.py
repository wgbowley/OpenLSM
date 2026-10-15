"""
Filename: optimizer.py

Description:
    converging random search using 1D analytical 
    simulation to refine slot size and slot thickness.
"""

import csv
from pathlib import Path
from random import random

from picounits import Q, Parser, DynamicLoader
from picounits import LENGTH, FORCE, CURRENT, POWER, NULLSET, TIME
from model.solver import Solver


# Loads unit system, material library & parameters
ROOT_DIR = Path(__file__).resolve().parents[0]

# Materials & Parameter files
parameters_path = ROOT_DIR / "parameters.uiv"
parameters = Parser.open(parameters_path, ROOT_DIR / "derived.ut")


def evaluate(param: DynamicLoader) -> tuple[Q, Q, Q]:
    """ Evaluates the proposed motor design over one pole pitch. """
    # Magnetic Solver for tubular linear motor
    solver = Solver(param)

    z_sample = param.numerics.displacement
    step_size = param.numerics.samples.displacement_size

    # Calculates the sampling domain & step sizes
    steps = int(z_sample / step_size)
    offset = - z_sample / 2

    # List to store
    force_result = []
    for index in range(0, steps):
        # Calculates the armature position
        z_output = index * step_size + offset
        z_solver = z_output - solver.armature_offset * LENGTH

        # Updates the currents & calculates force
        force = solver.compute_force(z_solver)

        force_result.append(force)

    # Computes the average force & force per ampere
    force_average = sum(force_result, 0 * FORCE) / len(force_result)
    force_per_amp = force_average / (solver.phase_rms_current * CURRENT)

    # Computes the force slew rate
    current_raise = solver.line_voltage / solver.phase_inductance
    force_per_second = force_per_amp * (current_raise * (CURRENT / TIME))

    # Returns the force per ampere and copper losses (first order losses)
    return force_average, force_per_amp, force_per_second, solver.average_losses * POWER


# Creates a optimization trace file
csv_filename = ROOT_DIR / "optimization_trace.csv"
csv_file = open(csv_filename, mode='w', newline='')
csv_writer = csv.writer(csv_file)

# Adds the specific rows the solver uses.
csv_writer.writerow([
    "Iteration", 
    "Axial Length (m)", 
    "Radial Thickness (m)",
    "Absolute Force (N)",
    "Force per Ampere (N/A)",
    "Force slew rate (N/s)"
    "Copper Loss (W)", 
    "Normalized Score"
])

# Optimization parameters
stall = 0
step_size = parameters.optimizer.size.start

# Series of best outputs
best_normalized        = - 1e10 * NULLSET
best_force             = - 1e10 * NULLSET
best_force_per_amp     = - 1e10 * FORCE / CURRENT
best_force_per_second  = - 1e10 * FORCE / TIME
best_losses            = - 1e10 * POWER
best_slot_axial_length = parameters.armature.slots.axial_length
best_slot_thickness    = parameters.armature.slots.radial_thickness

# Optimization loop
index = 0
for index in range(0, parameters.optimizer.samples.stripped):
    print(index, best_force, best_force_per_amp, best_force_per_second, best_losses)

    # Perturb current best geometry
    axial = abs(best_slot_axial_length + step_size * (random() * 2 - 1))
    thickness = abs(best_slot_thickness + step_size * (random() * 2 - 1))

    if axial > parameters.armature.slots.axial_pitch:
        # Skips invalid axial lengths
        continue

    # Sets new parameters to attempt and solves that set
    parameters.armature.slots.axial_length = axial
    parameters.armature.slots.radial_thickness = thickness

    try:
        average_force, force_per_amp, force_per_second, losses = evaluate(parameters)
    except:
        # Skips invalid designs
        continue

    # Normalizes the evaluation function
    normalized_F = average_force / parameters.optimizer.normalization.average_force
    normalized_F_A = force_per_amp / parameters.optimizer.normalization.force_per_ampere
    normalized_F_S = force_per_second / parameters.optimizer.normalization.force_per_second
    normalized_losses = losses / parameters.optimizer.normalization.copper_loss
    normalized = normalized_F + normalized_F_S + normalized_F_A - 3 * normalized_losses

    # Logs all evaluations
    csv_writer.writerow([
        index,
        axial.stripped,
        thickness.stripped,
        average_force.stripped,
        force_per_amp.stripped,
        force_per_second.stripped,
        losses.stripped,
        normalized.stripped
    ])

    if normalized > best_normalized:
        # Updates the best normalized & drops the stall
        best_normalized = normalized
        stall = 0

        # Saves the best outputs & inputs
        best_force = average_force
        best_force_per_amp = force_per_amp
        best_force_per_second = force_per_second
        best_losses = losses

        best_slot_axial_length = axial
        best_slot_thickness = thickness
        continue

    if stall > parameters.optimizer.stall_max.stripped:
        step_size /=2
        stall = 0
        print(f"Step size reduced to {step_size:.3f} due to stagnation.")

    if step_size < parameters.optimizer.size.min:
        print("Minimum step size reached. Ending optimization.")
        break

    stall += 1

csv_file.close()
print(f"Optimization script completed. stall = {stall}, size={step_size}, index={index}")
