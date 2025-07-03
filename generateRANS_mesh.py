import gmsh

gmsh.initialize()
gmsh.model.add("wing_mesh")

# Wing geometry (m)
c, b, t = 0.04, 0.40427, 0.002

# Wind-tunnel extents (5 c upstream, 7.5 c downstream, 0.5 b margin)
min_x, max_x = -0.2, 0.3
min_y, max_y = -b/2 - 0.2, b/2 + 0.2
min_z, max_z = -0.2, 0.2

# 1) Create domain & wing
domain = gmsh.model.occ.addBox(min_x, min_y, min_z,
                               max_x - min_x,
                               max_y - min_y,
                               max_z - min_z)
wing   = gmsh.model.occ.addBox(0.0, -b/2, -t/2,
                               c, b, t)

# 2) Carve out the wing
gmsh.model.occ.cut([(3, domain)], [(3, wing)],
                   removeObject=True,
                   removeTool=True)
gmsh.model.occ.synchronize()

# 3) Grab the only remaining volume = fluid
fluid_tag = gmsh.model.getEntities(dim=3)[0][1]

# 4) Identify surfaces for physical groups
inlet_surf     = []
outlet_surf    = []
horizontal_surf = []
airfoil_surf   = []

for _, tag in gmsh.model.getEntities(dim=2):
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)

    # inlet / outlet
    if abs(xmin - min_x) < 1e-6 and abs(xmax - min_x) < 1e-6:
        inlet_surf.append(tag)
    elif abs(xmin - max_x) < 1e-6 and abs(xmax - max_x) < 1e-6:
        outlet_surf.append(tag)

    # walls (all four) → horizontal
    elif (abs(zmin - min_z) < 1e-6 and abs(zmax - min_z) < 1e-6) or \
         (abs(zmin - max_z) < 1e-6 and abs(zmax - max_z) < 1e-6) or \
         (abs(ymin - min_y) < 1e-6 and abs(ymax - min_y) < 1e-6) or \
         (abs(ymin - max_y) < 1e-6 and abs(ymax - max_y) < 1e-6):
        horizontal_surf.append(tag)

    # the six wing faces → airfoil
    else:
        airfoil_surf.append(tag)

# 5) Create your five physical groups
gmsh.model.addPhysicalGroup(2, inlet_surf,     1)
gmsh.model.setPhysicalName(2, 1, "inlet")

gmsh.model.addPhysicalGroup(2, outlet_surf,    2)
gmsh.model.setPhysicalName(2, 2, "outlet")

gmsh.model.addPhysicalGroup(2, horizontal_surf,3)
gmsh.model.setPhysicalName(2, 3, "horizontal")

gmsh.model.addPhysicalGroup(2, airfoil_surf,   4)
gmsh.model.setPhysicalName(2, 4, "airfoil")

gmsh.model.addPhysicalGroup(3, [fluid_tag],    5)
gmsh.model.setPhysicalName(3, 5, "fluid")

 # ─────────────────────────────────────────────────────────────────────────────
 # 6) Refined multi-stage grading via TWO Thresholds + Min field
# Stage sizes and distances (more aggressive)
size_near, size_mid, size_far = 0.001, 0.01, 0.01    # 1 mm → 100 mm → 200 mm
dist_min,  dist_mid,  dist_max = 0.001, 0.01,  0.01   # 3 mm → 150 mm → 500 mm

 # 6a) Distance field from the wing surfaces
fldD = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(fldD, "SurfacesList", airfoil_surf)
gmsh.model.mesh.field.setNumber(fldD, "Sampling", 50)

# 6b) First threshold: 1 mm → size_mid over dist_min → dist_mid
t1 = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(t1, "IField",   fldD)
gmsh.model.mesh.field.setNumber(t1, "LcMin",    size_near)
gmsh.model.mesh.field.setNumber(t1, "LcMax",    size_mid)
gmsh.model.mesh.field.setNumber(t1, "DistMin",  dist_min)
gmsh.model.mesh.field.setNumber(t1, "DistMax",  dist_mid)

# 6c) Second threshold: size_mid → size_far over dist_mid → dist_max
t2 = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(t2, "IField",   fldD)
gmsh.model.mesh.field.setNumber(t2, "LcMin",    size_mid)
gmsh.model.mesh.field.setNumber(t2, "LcMax",    size_far)
gmsh.model.mesh.field.setNumber(t2, "DistMin",  dist_mid)
gmsh.model.mesh.field.setNumber(t2, "DistMax",  dist_max)

# 6d) Combine them: pick the minimum size at each point
fmin = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", [t1, t2])
gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
# ─────────────────────────────────────────────────────────────────────────────


# 7) Enable a little optimization and size bounds
gmsh.option.setNumber("Mesh.Optimize", 1)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", size_near)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", size_far)

# 8) Generate & write
gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("wing_mesh.msh")
gmsh.finalize()
