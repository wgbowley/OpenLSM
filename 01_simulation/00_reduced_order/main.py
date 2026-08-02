import numpy as np
import matplotlib.pyplot as plt
from model.physics.field import compute_dipole_field_strength

# Parameters
translation = 0.05
current = 1.0  # Amperes
turns = 100    # Number of turns
length = 0.05  # meters (50mm)
radius = 0.01  # meters (10mm)

# Generate z positions from -0.2 to 0.2
z_positions = np.linspace(-0.2, 0.2, 1000)

# Compute field strength at each position
field_strengths = [compute_dipole_field_strength(z, 0.1, 10, 0.01)
                   for z in z_positions]

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(z_positions, field_strengths, linewidth=2)
plt.xlabel('Z Position (m)', fontsize=12)
plt.ylabel('Field Strength (A/m)', fontsize=12)
plt.title('Magnetic Field Strength Along Z-Axis', fontsize=14)
plt.grid(True, alpha=0.3)
plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
plt.xlim(-0.2, 0.2)
plt.tight_layout()
plt.show()