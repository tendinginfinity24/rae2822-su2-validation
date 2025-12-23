import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np

# --- 1. NASA Case 6 Experimental Data (Cook et al.) ---
exp_x_c = [
    0.006, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095,
    0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.175, 0.195, 0.215, 0.235,
    0.255, 0.275, 0.295, 0.315, 0.335, 0.355, 0.375, 0.395, 0.415, 0.435,
    0.455, 0.475, 0.495, 0.515, 0.535, 0.555, 0.575, 0.595, 0.615, 0.635,
    0.655, 0.675, 0.695, 0.715, 0.735, 0.755, 0.775, 0.795, 0.815, 0.835,
    0.855, 0.875, 0.895, 0.915, 0.935, 0.950
]

exp_cp = [
    -0.457, -0.650, -0.748, -0.806, -0.849, -0.887, -0.916, -0.942, -0.963, -0.985,
    -0.999, -1.015, -1.026, -1.037, -1.047, -1.056, -1.077, -1.096, -1.118, -1.139,
    -1.161, -1.185, -1.213, -1.241, -1.272, -1.298, -1.275, -1.233, -1.203, -1.184,
    -1.166, -1.157, -1.159, -1.171, -1.192, -1.216, -0.677, -0.584, -0.528, -0.479,
    -0.428, -0.372, -0.320, -0.271, -0.222, -0.177, -0.133, -0.091, -0.046, -0.003,
    0.040, 0.080, 0.120, 0.160, 0.196, 0.220
]

# --- 2. Load SU2 Results with PyVista ---
filename = "surface_flow.vtu"

try:
    mesh = pv.read(filename)
except FileNotFoundError:
    print(f"Error: Could not find '{filename}'. Did you run SU2 yet?")
    exit()

# Extract Data Arrays
x = mesh.points[:, 0]  # X coordinates
y = mesh.points[:, 1]  # Y coordinates (needed for sorting upper/lower)
cp = mesh.point_data["Pressure_Coefficient"]

# --- 3. Sort Data for Plotting ---
# Since VTU is unstructured, we separate Upper/Lower surfaces to get clean lines
# Assuming Airfoil Chord is along X-axis
mask_upper = y >= 0
mask_lower = y < 0

# Create structured arrays for plotting
upper_x = x[mask_upper]
upper_cp = cp[mask_upper]
lower_x = x[mask_lower]
lower_cp = cp[mask_lower]

# Sort by x coordinate so the line plot connects dots in order
# (argsort returns the indices that would sort the array)
sort_idx_upper = np.argsort(upper_x)
sort_idx_lower = np.argsort(lower_x)

# Apply sorting
upper_x = upper_x[sort_idx_upper]
upper_cp = upper_cp[sort_idx_upper]
lower_x = lower_x[sort_idx_lower]
lower_cp = lower_cp[sort_idx_lower]

# --- 4. Plotting ---
plt.figure(figsize=(10, 6))

# Plot SU2 Data (Now as two clean lines)
# We plot Upper and Lower separately so the line doesn't wrap around the nose/tail incorrectly
plt.plot(upper_x, upper_cp, label="SU2 Upper", color="red", linewidth=2)
plt.plot(lower_x, lower_cp, label="SU2 Lower", color="blue", linewidth=2, linestyle="--")

# Plot Experimental Data
plt.scatter(exp_x_c, exp_cp, label="NASA Case 6 (Exp)", color="black", marker="o", s=40, zorder=5)

# Formatting
plt.gca().invert_yaxis()  # Negative Cp up
plt.title(f"RAE 2822 Case 6 Validation\nMa=0.729, Re=6.5E6", fontsize=14)
plt.xlabel("x/c", fontsize=12)
plt.ylabel("$C_p$", fontsize=12)
plt.grid(True, linestyle="--", alpha=0.7)
plt.legend()

# Save and Show
plt.savefig("cp_validation_vtu.png", dpi=300)
plt.show()