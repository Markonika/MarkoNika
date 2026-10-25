#!/usr/bin/env python3
"""
Thin 3D wedge mesh of the EVOL/MSFR curved ("torus") core cavity -- v2.

Corrected per user feedback on v1 (generate_mesh.py):
  - The core body is now a single smooth, closed, lens/"stadium"-shaped
    cross-section (like a watermelon-wedge rind) -- no flat rectangular
    duct boxes bolted onto it at right angles.
  - Inlet and outlet are genuine circular pipes (cylinders), fused onto
    the curved body with a 3D boolean union, not rectangular slots.
  - There is exactly one axis of revolution (the true vertical reactor
    axis, r=0), and exactly two flat "symmetry" cut faces bounding the
    thin wedge -- necessary for a thin-wedge axisymmetric-equivalent
    model at all (a wedge cut from a round body always has two flat
    cut faces unless the wedge spans the full 360 degrees), but no
    rectangular geometry anywhere else.
  - The wedge angle is no longer an arbitrary thin 4 degrees: a real
    circular pipe doesn't fit inside an arbitrarily thin angular slice,
    so the angle is computed here to just contain the pipe's diameter
    (with margin) at the radius where it attaches.

SOURCED vs NOT SOURCED (unchanged from v1 except as noted):
  - Core envelope (R_MIN, R_MAX, H_IN, H_OUT): Fang et al. (2024),
    Nucl. Eng. Des. 424, 113294 -- Table 1 / Section 2.
  - Profile curve shape: engineering approximation within that sourced
    envelope (smooth cosine-eased curve), not a pixel-trace of the
    published CAD/mesh.
  - PIPE_DIAMETER, PIPE_LENGTH: NOT sourced from any paper reviewed.
    Flagged placeholders (same order of magnitude as the duct height
    used in v1, now reinterpreted as a true pipe diameter).
"""
import gmsh
import numpy as np

# ---------------------------------------------------------------------
# Sourced envelope (Fang et al. 2024)
# ---------------------------------------------------------------------
R_MIN = 1.054
R_MAX = 1.53
H_IN = 1.6
H_OUT = 2.65

half_in = H_IN / 2
half_out = H_OUT / 2

# ---------------------------------------------------------------------
# NOT sourced -- placeholder pipe geometry, flagged explicitly
# ---------------------------------------------------------------------
PIPE_DIAMETER = 0.30   # m
PIPE_LENGTH = 0.40     # m, how far the pipe stub extends beyond R_MAX
Z_INLET_CENTER = -1.0  # m, axial position of inlet pipe centerline
Z_OUTLET_CENTER = 1.0  # m, axial position of outlet pipe centerline
PIPE_MARGIN = 1.5      # safety margin on the wedge angle vs. bare pipe fit

MESH_SIZE_WALL = 0.05
MESH_SIZE_PIPE = 0.03

# wedge angle: must be wide enough that the pipe (diameter PIPE_DIAMETER,
# attached at radius ~R_MAX) doesn't poke through the two angular cut
# faces. At radius r the arc width per unit angle is r*theta, so require
# R_MAX * wedge_angle >= PIPE_DIAMETER * PIPE_MARGIN.
wedge_angle = PIPE_MARGIN * PIPE_DIAMETER / R_MAX
print(f"Computed wedge angle: {np.degrees(wedge_angle):.2f} deg "
      f"(vs. an arbitrary 4 deg in v1 -- now sized to fit the pipe)")
theta_mid = wedge_angle / 2


def smooth_T(r):
    s = (r - R_MIN) / (R_MAX - R_MIN)
    return half_in + (half_out - half_in) * 0.5 * (1 - np.cos(np.pi * s))


gmsh.initialize()
gmsh.model.add("msfr_wedge_v2")
occ = gmsh.model.occ

# ===================== CORE BODY: smooth closed profile, revolved =====
# Profile in the (x=axial, y=radial) plane at the wedge's theta=0 face.
N = 16
pts = []
pts.append(occ.addPoint(-half_in, R_MIN, 0.0, MESH_SIZE_WALL))  # inner-bottom
pts.append(occ.addPoint(half_in, R_MIN, 0.0, MESH_SIZE_WALL))   # inner-top
for s in np.linspace(0, 1, N)[1:-1]:
    r = R_MIN + s * (R_MAX - R_MIN)
    pts.append(occ.addPoint(smooth_T(r), r, 0.0, MESH_SIZE_WALL))
pts.append(occ.addPoint(half_out, R_MAX, 0.0, MESH_SIZE_WALL))  # outer-top
pts.append(occ.addPoint(-half_out, R_MAX, 0.0, MESH_SIZE_WALL))  # outer-bottom
for s in np.linspace(0, 1, N)[1:-1][::-1]:
    r = R_MIN + s * (R_MAX - R_MIN)
    pts.append(occ.addPoint(-smooth_T(r), r, 0.0, MESH_SIZE_WALL))

lines = [occ.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
loop = occ.addCurveLoop(lines)
profile_surf = occ.addPlaneSurface([loop])
occ.synchronize()

core = occ.revolve([(2, profile_surf)], 0, 0, 0, 1, 0, 0, wedge_angle)
core_vol = [e for e in core if e[0] == 3][0]

# ===================== PIPES: real circular cylinders =====================
def add_pipe(z_center):
    # base point slightly inside the core body (so the boolean union
    # has solid overlap, not a tangent/zero-thickness join), extending
    # outward to R_MAX + PIPE_LENGTH, oriented along the wedge's
    # mid-angle radial direction.
    r0 = R_MIN  # start from deep inside the body for a robust union
    length = (R_MAX + PIPE_LENGTH) - r0
    base = (z_center, r0 * np.cos(theta_mid), r0 * np.sin(theta_mid))
    direction = (0.0, length * np.cos(theta_mid), length * np.sin(theta_mid))
    tag = occ.addCylinder(*base, *direction, PIPE_DIAMETER / 2)
    return tag


inlet_pipe = add_pipe(Z_INLET_CENTER)
outlet_pipe = add_pipe(Z_OUTLET_CENTER)
occ.synchronize()

fused, _ = occ.fuse([core_vol], [(3, inlet_pipe), (3, outlet_pipe)])
occ.synchronize()

fused_vol = fused[0][1]

# ===================== CLASSIFY BOUNDARY FACES =====================
boundary = gmsh.model.getBoundary([(3, fused_vol)], oriented=False, recursive=False)
surf_tags = [t for d, t in boundary]


def face_com(tag):
    return np.array(occ.getCenterOfMass(2, tag))


def pipe_end_target(z_center):
    r_end = R_MAX + PIPE_LENGTH
    return np.array([z_center, r_end * np.cos(theta_mid), r_end * np.sin(theta_mid)])


inlet_target = pipe_end_target(Z_INLET_CENTER)
outlet_target = pipe_end_target(Z_OUTLET_CENTER)

sym0_normal = np.array([0.0, 0.0, 1.0])     # theta=0 face: z_gmsh = 0 plane
sym1_normal = np.array([0.0, -np.sin(wedge_angle), np.cos(wedge_angle)])

# A true pipe end cap is a flat disk of area pi*(D/2)^2 -- the
# lateral (cylindrical) surface of the pipe stub has a much larger
# area and must NOT be classified as inlet/outlet just because its
# centroid is nearby.
expected_cap_area = np.pi * (PIPE_DIAMETER / 2) ** 2

classified = {}
remaining = []
for tag in surf_tags:
    com = face_com(tag)
    area = occ.getMass(2, tag)
    is_cap_area = abs(area - expected_cap_area) < 0.1 * expected_cap_area
    if is_cap_area and np.linalg.norm(com - inlet_target) < PIPE_DIAMETER:
        classified[tag] = "inlet"
    elif is_cap_area and np.linalg.norm(com - outlet_target) < PIPE_DIAMETER:
        classified[tag] = "outlet"
    else:
        remaining.append(tag)

# symmetry faces: planar, with z_gmsh (for sym0) or the sym1-plane
# coordinate ~ 0 for essentially the whole face. Identify by checking
# that every vertex of the face lies (numerically) on the plane.
def on_plane(tag, point_on_plane, normal, tol=1e-4):
    pts_on_face = gmsh.model.getBoundary([(2, tag)], oriented=False, recursive=True)
    for d, p in pts_on_face:
        xyz = np.array(gmsh.model.getValue(0, p, []))
        if abs(np.dot(xyz - point_on_plane, normal)) > tol:
            return False
    return True


wall_tags = []
for tag in remaining:
    if on_plane(tag, np.array([0.0, 0.0, 0.0]), sym0_normal):
        classified[tag] = "symmetry_0"
    elif on_plane(tag, np.array([0.0, 0.0, 0.0]), sym1_normal):
        classified[tag] = "symmetry_1"
    else:
        wall_tags.append(tag)

inlet_faces = [t for t, n in classified.items() if n == "inlet"]
outlet_faces = [t for t, n in classified.items() if n == "outlet"]
sym0_faces = [t for t, n in classified.items() if n == "symmetry_0"]
sym1_faces = [t for t, n in classified.items() if n == "symmetry_1"]

print(f"inlet faces: {inlet_faces}, outlet faces: {outlet_faces}")
print(f"symmetry_0 faces: {sym0_faces}, symmetry_1 faces: {sym1_faces}")
print(f"wall faces: {len(wall_tags)}")

gmsh.model.addPhysicalGroup(3, [fused_vol], name="fluid")
gmsh.model.addPhysicalGroup(2, inlet_faces, name="inlet")
gmsh.model.addPhysicalGroup(2, outlet_faces, name="outlet")
gmsh.model.addPhysicalGroup(2, sym0_faces, name="symmetry_0")
gmsh.model.addPhysicalGroup(2, sym1_faces, name="symmetry_1")
gmsh.model.addPhysicalGroup(2, wall_tags, name="wall")

# finer mesh on the pipes
gmsh.model.mesh.setSize(
    gmsh.model.getBoundary([(2, t) for t in inlet_faces + outlet_faces],
                            recursive=True),
    MESH_SIZE_PIPE)

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 8)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)

gmsh.model.mesh.generate(3)

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("/home/user/MarkoNika/msfr_wedge/msfr_wedge_v2.msh")

n_nodes = len(gmsh.model.mesh.getNodes()[0])
n_elems = sum(len(e) for e in gmsh.model.mesh.getElements(3)[1])
print(f"Mesh written: {n_nodes} nodes, {n_elems} tetrahedra")

gmsh.finalize()
