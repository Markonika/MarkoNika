#!/usr/bin/env python3
"""
Thin 3D wedge mesh of the EVOL/MSFR curved ("torus") core cavity, for
axisymmetric-equivalent CFD in Fluidity.

Fluidity has no native axisymmetric (r-z) solver, so a true 2D-axisymmetric
run (as used in Fluent/Nek5000 in the literature) isn't available. Instead
this builds a thin azimuthal wedge of the real 3D geometry -- exactly the
technique used by Fang et al. (2024) ("3-D wedge domains") and Laureau
et al. (2022) ("1/32nd MSFR model") -- with symmetry (no_normal_flow) BCs
on the two angular cut faces. This is ordinary 3D Navier-Stokes with no
special solver support needed, and is physically correct in the limit of
a small wedge angle (here 4 degrees).

GEOMETRY SOURCE AND HONESTY NOTE
---------------------------------
The overall envelope is sourced directly from Fang, Tano, Saini, Tomboulides,
Coppo-Leite, Merzari, Feng, Shaver (2024), "CFD simulations of Molten Salt
Fast Reactor core cavity flows," Nuclear Engineering and Design 424, 113294:
    - core radius ranges 1.05 - 1.53 m (Table 1: minimum core radius 1.054 m)
    - core length 1.6 m along the centerline, 2.65 m at the periphery
This defines the (r, z) bounding envelope of the curved "torus" cavity
("Geometry II" of Rouch et al. 2014, which Fang et al. adopt as their
reference geometry). Neither paper publishes the exact analytic/CAD curve
of the cavity wall. The profile built here is therefore an ENGINEERING
APPROXIMATION within that sourced envelope -- a smooth cosine-eased curve
from the narrow inner cap (height 1.6 m) to the taller outer rim (height
2.65 m) -- not a pixel-trace of the published CAD/mesh figures. The inlet
and outlet duct dimensions (0.3 m tall, protruding 0.4 m) are NOT sourced
from any paper and are a reasonable placeholder, flagged here explicitly;
change DUCT_HEIGHT / DUCT_LENGTH below if better data becomes available.

Everything else -- flow direction (inlet at the bottom of the outer wall,
outlet at the top of the outer wall), operating conditions, and fluid
properties -- is sourced; see the accompanying case README.
"""
import gmsh
import numpy as np

# ---------------------------------------------------------------------
# Sourced envelope (Fang et al. 2024)
# ---------------------------------------------------------------------
R_MIN = 1.054       # m, minimum core radius (Table 1)
R_MAX = 1.53        # m, maximum core radius (Section 2, rounded)
H_IN = 1.6          # m, core length along the centerline (-> half = 0.8)
H_OUT = 2.65        # m, core length at the periphery (-> half = 1.325)

# ---------------------------------------------------------------------
# NOT sourced -- placeholder duct geometry, flagged explicitly
# ---------------------------------------------------------------------
DUCT_HEIGHT = 0.30  # m, inlet/outlet duct opening height
DUCT_LENGTH = 0.40  # m, radial protrusion of the duct beyond R_MAX

WEDGE_ANGLE_DEG = 4.0   # thin azimuthal wedge (axisymmetric-equivalent)
N_PROFILE_PTS = 14      # points along each curved top/bottom wall segment
MESH_SIZE_WALL = 0.05   # m, target element size near walls
MESH_SIZE_BULK = 0.10   # m, target element size in the bulk

half_in = H_IN / 2
half_out = H_OUT / 2
z_duct_top_inner = half_out - DUCT_HEIGHT   # inner (smaller |z|) edge of the duct opening
z_duct_bot_inner = -z_duct_top_inner

gmsh.initialize()
gmsh.model.add("msfr_wedge")
geo = gmsh.model.geo


def smooth_T(r):
    """Cosine-eased half-height of the core cavity at radius r."""
    s = (r - R_MIN) / (R_MAX - R_MIN)
    return half_in + (half_out - half_in) * 0.5 * (1 - np.cos(np.pi * s))


# ---------------------------------------------------------------------
# Build the closed (r, z) profile, going clockwise from the bottom of
# the inner cap, with z plotted as gmsh's y-coordinate (profile lives
# in the x-y plane at z_gmsh=0; revolve is about the x-axis, which
# plays the role of the true vertical reactor axis).
# ---------------------------------------------------------------------
pts = []


def add_pt(r, z, mesh_size):
    # Embed (r, z) as gmsh (x, y) = (z, r), so that revolving about the
    # gmsh x-axis below correctly keeps the physical axial coordinate z
    # fixed while sweeping the physical radius r around it -- i.e. the
    # gmsh x-axis IS the true vertical reactor axis. (Embedding as
    # (x, y) = (r, z) and revolving about x would instead keep r fixed
    # and sweep z through the origin, self-intersecting the geometry.)
    tag = geo.addPoint(z, r, 0.0, mesh_size)
    pts.append(tag)
    return tag

# 1. bottom of inner cap
p_inner_bot = add_pt(R_MIN, -half_in, MESH_SIZE_WALL)
# 2. top of inner cap
p_inner_top = add_pt(R_MIN, half_in, MESH_SIZE_WALL)

# top curved wall: r from R_MIN to R_MAX, z = +T(r)
top_curve_pts = [p_inner_top]
for s in np.linspace(0, 1, N_PROFILE_PTS)[1:-1]:
    r = R_MIN + s * (R_MAX - R_MIN)
    top_curve_pts.append(add_pt(r, smooth_T(r), MESH_SIZE_WALL))
p_top_outer = add_pt(R_MAX, half_out, MESH_SIZE_WALL)
top_curve_pts.append(p_top_outer)

# outlet duct (top): out, down (outlet face), in
p_outlet_out_top = add_pt(R_MAX + DUCT_LENGTH, half_out, MESH_SIZE_WALL)
p_outlet_out_bot = add_pt(R_MAX + DUCT_LENGTH, z_duct_top_inner, MESH_SIZE_WALL)
p_outlet_in = add_pt(R_MAX, z_duct_top_inner, MESH_SIZE_WALL)

# main outer (equatorial) wall: down from outlet duct to inlet duct
p_inlet_in = add_pt(R_MAX, z_duct_bot_inner, MESH_SIZE_WALL)

# inlet duct (bottom): out, down (inlet face), in
p_inlet_out_top = add_pt(R_MAX + DUCT_LENGTH, z_duct_bot_inner, MESH_SIZE_WALL)
p_inlet_out_bot = add_pt(R_MAX + DUCT_LENGTH, -half_out, MESH_SIZE_WALL)
p_bottom_outer = add_pt(R_MAX, -half_out, MESH_SIZE_WALL)

# bottom curved wall: r from R_MAX back to R_MIN, z = -T(r)
bot_curve_pts = [p_bottom_outer]
for s in np.linspace(0, 1, N_PROFILE_PTS)[1:-1][::-1]:
    r = R_MIN + s * (R_MAX - R_MIN)
    bot_curve_pts.append(add_pt(r, -smooth_T(r), MESH_SIZE_WALL))
bot_curve_pts.append(p_inner_bot)

geo.synchronize()

# ---------------------------------------------------------------------
# Curves
# ---------------------------------------------------------------------
def poly_lines(point_list):
    """Approximate a smooth curve with straight segments between
    consecutive points (avoids NURBS ruled-surface degeneracies when
    later revolved through a thin wedge angle)."""
    return [geo.addLine(point_list[i], point_list[i + 1])
            for i in range(len(point_list) - 1)]


l_inner_cap = geo.addLine(p_inner_bot, p_inner_top)
l_top_wall = poly_lines(top_curve_pts)
l_outlet_top = geo.addLine(p_top_outer, p_outlet_out_top)
l_outlet_face = geo.addLine(p_outlet_out_top, p_outlet_out_bot)
l_outlet_bot = geo.addLine(p_outlet_out_bot, p_outlet_in)
l_equator = geo.addLine(p_outlet_in, p_inlet_in)
l_inlet_top = geo.addLine(p_inlet_in, p_inlet_out_top)
l_inlet_face = geo.addLine(p_inlet_out_top, p_inlet_out_bot)
l_inlet_bot = geo.addLine(p_inlet_out_bot, p_bottom_outer)
l_bot_wall = poly_lines(bot_curve_pts)

loop = geo.addCurveLoop([
    l_inner_cap, *l_top_wall, l_outlet_top, l_outlet_face, l_outlet_bot,
    l_equator, l_inlet_top, l_inlet_face, l_inlet_bot, *l_bot_wall,
])
surf = geo.addPlaneSurface([loop])
geo.synchronize()

# ---------------------------------------------------------------------
# Revolve about the x-axis (r=0, the true vertical reactor axis) by a
# thin wedge angle to get the axisymmetric-equivalent 3D domain.
# ---------------------------------------------------------------------
angle = np.deg2rad(WEDGE_ANGLE_DEG)
out = geo.revolve([(2, surf)], 0, 0, 0, 1, 0, 0, angle)
geo.synchronize()

# out is a list of (dim, tag) entries: [volume, top-face(symmetry_1),
# then one entry per swept boundary curve of the original profile, in
# the same order as the curve loop above.
vol_tag = [e[1] for e in out if e[0] == 3][0]
swept_faces = [e[1] for e in out if e[0] == 2]
# the first swept face corresponds to the "end cap" (symmetry_1); the
# rest correspond, in order, to the profile curves above
symmetry_1_tag = swept_faces[0]
swept_wall_faces = swept_faces[1:]

curve_order = [l_inner_cap, *l_top_wall, l_outlet_top, l_outlet_face,
               l_outlet_bot, l_equator, l_inlet_top, l_inlet_face,
               l_inlet_bot, *l_bot_wall]
face_by_curve = dict(zip(curve_order, swept_wall_faces))

outlet_face = face_by_curve[l_outlet_face]
inlet_face = face_by_curve[l_inlet_face]
wall_faces = [face_by_curve[c] for c in curve_order
              if c not in (l_outlet_face, l_inlet_face)]

gmsh.model.addPhysicalGroup(3, [vol_tag], name="fluid")
gmsh.model.addPhysicalGroup(2, [surf], name="symmetry_0")
gmsh.model.addPhysicalGroup(2, [symmetry_1_tag], name="symmetry_1")
gmsh.model.addPhysicalGroup(2, [inlet_face], name="inlet")
gmsh.model.addPhysicalGroup(2, [outlet_face], name="outlet")
gmsh.model.addPhysicalGroup(2, wall_faces, name="wall")

# ---------------------------------------------------------------------
# Mesh size: use the explicit per-point characteristic lengths set
# above (uniform-ish near walls, coarser in the bulk via the few bulk
# points); a Distance/Threshold background field combined with the
# thin wedge angle produced degenerate near-duct triangles, so it's
# left out here in favour of the simpler, more robust point sizing.
# ---------------------------------------------------------------------
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay

gmsh.model.mesh.generate(3)

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("/home/user/MarkoNika/msfr_wedge/msfr_wedge.msh")

n_nodes = len(gmsh.model.mesh.getNodes()[0])
n_elems = sum(len(e) for e in gmsh.model.mesh.getElements(3)[1])
print(f"Mesh written: {n_nodes} nodes, {n_elems} tetrahedra")

gmsh.finalize()
