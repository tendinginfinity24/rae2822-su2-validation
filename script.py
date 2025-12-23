import gmsh
import sys

# --- 1. Parameters (RAE 2822 Case 6 Calibration) ---
CHORD = 1.0
lc_farfield = 5.0
lc_airfoil = 0.003  # Finer general distribution
lc_le_te = 0.0008   # Sharp LE and TE resolution
lc_shock = 0.001    # Refinement for x/c 0.4 to 0.7

# Boundary Layer (Targeting y+ ~ 1)
first_layer = 2.5e-6 
growth_rate = 1.18
num_layers = 35

# --- 2. Initialize ---
gmsh.initialize()
gmsh.model.add("RAE2822_Case6_Validation")

# --- 3. Airfoil Coordinates ---
rae_upper = [
    (0.00000, 0.00000), (0.00500, 0.01085), (0.01250, 0.01569), (0.02500, 0.02081),
    (0.05000, 0.02856), (0.07500, 0.03483), (0.10000, 0.04036), (0.15000, 0.04987),
    (0.20000, 0.05770), (0.25000, 0.06394), (0.30000, 0.06861), (0.35000, 0.07174),
    (0.40000, 0.07328), (0.45000, 0.07328), (0.50000, 0.07179), (0.55000, 0.06893),
    (0.60000, 0.06484), (0.65000, 0.05961), (0.70000, 0.05331), (0.75000, 0.04603),
    (0.80000, 0.03787), (0.85000, 0.02897), (0.90000, 0.01948), (0.95000, 0.00960),
    (0.97500, 0.00465), (1.00000, 0.00000)
]
rae_lower = [
    (0.00000, 0.00000), (0.00500, -0.00105), (0.01250, -0.00213), (0.02500, -0.00392),
    (0.05000, -0.00742), (0.07500, -0.01083), (0.10000, -0.01416), (0.15000, -0.02069),
    (0.20000, -0.02693), (0.25000, -0.03260), (0.30000, -0.03758), (0.35000, -0.04173),
    (0.40000, -0.04491), (0.45000, -0.04706), (0.50000, -0.04812), (0.55000, -0.04812),
    (0.60000, -0.04714), (0.65000, -0.04533), (0.70000, -0.04273), (0.75000, -0.03930),
    (0.80000, -0.03493), (0.85000, -0.02941), (0.90000, -0.02242), (0.95000, -0.01342),
    (0.97500, -0.00766), (1.00000, 0.00000)
]

# --- 4. Geometry Construction ---
p_le = gmsh.model.geo.addPoint(0, 0, 0, lc_le_te)
p_te = gmsh.model.geo.addPoint(1, 0, 0, lc_le_te)

pts_u = [p_le]
for (x, y) in rae_upper[1:-1]:
    # Refine specifically in the shock region (0.4 < x < 0.7)
    lc = lc_shock if 0.4 < x < 0.7 else lc_airfoil
    pts_u.append(gmsh.model.geo.addPoint(x, y, 0, lc))
pts_u.append(p_te)

pts_l = [p_te]
for (x, y) in reversed(rae_lower[1:-1]):
    pts_l.append(gmsh.model.geo.addPoint(x, y, 0, lc_airfoil))
pts_l.append(p_le)

l_upper = gmsh.model.geo.addSpline(pts_u)
l_lower = gmsh.model.geo.addSpline(pts_l)

# Farfield
radius = 30.0
c_x, c_y = 0.5, 0.0
p_c = gmsh.model.geo.addPoint(c_x, c_y, 0, lc_farfield)
pf1 = gmsh.model.geo.addPoint(c_x + radius, c_y, 0, lc_farfield)
pf2 = gmsh.model.geo.addPoint(c_x, c_y + radius, 0, lc_farfield)
pf3 = gmsh.model.geo.addPoint(c_x - radius, c_y, 0, lc_farfield)
pf4 = gmsh.model.geo.addPoint(c_x, c_y - radius, 0, lc_farfield)

arc1 = gmsh.model.geo.addCircleArc(pf1, p_c, pf2)
arc2 = gmsh.model.geo.addCircleArc(pf2, p_c, pf3)
arc3 = gmsh.model.geo.addCircleArc(pf3, p_c, pf4)
arc4 = gmsh.model.geo.addCircleArc(pf4, p_c, pf1)

loop_a = gmsh.model.geo.addCurveLoop([l_upper, l_lower])
loop_f = gmsh.model.geo.addCurveLoop([arc1, arc2, arc3, arc4])
surf = gmsh.model.geo.addPlaneSurface([loop_f, loop_a])

gmsh.model.geo.synchronize()

# --- 5. Boundary Layer Setup ---
f = gmsh.model.mesh.field.add("BoundaryLayer")
gmsh.model.mesh.field.setNumbers(f, "CurvesList", [l_upper, l_lower])
gmsh.model.mesh.field.setNumber(f, "Size", first_layer)
gmsh.model.mesh.field.setNumber(f, "Ratio", growth_rate)
gmsh.model.mesh.field.setNumber(f, "Thickness", 0.03) # Resolves the full BL
gmsh.model.mesh.field.setNumber(f, "Quads", 1)       # Crucial: generates Quads
gmsh.model.mesh.field.setAsBoundaryLayer(f)

# --- 6. Global Mesh Settings ---
gmsh.option.setNumber("Mesh.Algorithm", 6) # Frontal-Delaunay for 2D
gmsh.option.setNumber("Mesh.Smoothing", 10)

# --- 7. Finalize and Generate ---
gmsh.model.addPhysicalGroup(1, [l_upper, l_lower], name="AIRFOIL")
gmsh.model.addPhysicalGroup(1, [arc1, arc2, arc3, arc4], name="FARFIELD")
gmsh.model.addPhysicalGroup(2, [surf], name="FLUID")

# Define a box for wake refinement
f_wake = gmsh.model.mesh.field.add("Box")
gmsh.model.mesh.field.setNumber(f_wake, "VIn", 0.005)    # Size inside the box
gmsh.model.mesh.field.setNumber(f_wake, "VOut", 0.5)     # Size outside
gmsh.model.mesh.field.setNumber(f_wake, "XMin", 1.0)     # Starts at Trailing Edge
gmsh.model.mesh.field.setNumber(f_wake, "XMax", 4.0)     # Extends 3 chords back
gmsh.model.mesh.field.setNumber(f_wake, "YMin", -0.2)
gmsh.model.mesh.field.setNumber(f_wake, "YMax", 0.2)

# Use the 'Min' field to combine the Boundary Layer and the Wake
f_final = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(f_final, "FieldsList", [f, f_wake]) # 'f' is your BL field
gmsh.model.mesh.field.setAsBackgroundMesh(f_final)

gmsh.model.mesh.generate(2)
gmsh.write("rae2822_case6.su2")

if 'close' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()