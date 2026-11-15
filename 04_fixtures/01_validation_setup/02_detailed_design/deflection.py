""" 
Filename: Deflection.py 

Description:
    Calculates the deflection of two parallel 
    plates under cantilever bending load.
"""

from picounits import resolve_derived
from picounits import LENGTH, FORCE, PRESSURE, MEGA, MILLI

# Resolves derived units from path
resolve_derived()

# Parameters
youngs_modules = 5000 * MEGA * PRESSURE

height = 15 * MILLI * LENGTH
width = 45 * MILLI * LENGTH
thickness = 3 * MILLI * LENGTH

# Calculates the inertia of the cross section
inertia =  width * thickness ** 3 / 12

load = 1 * FORCE

for i in range(1, 20):
    # Assumes: 1/2 of the deflection of a single plate
    deflection = (load * height ** 3) / (3 * youngs_modules * inertia)
    deflection /= 2

    print(f"{load}: {deflection:.3f}")
    load += 1 * FORCE
