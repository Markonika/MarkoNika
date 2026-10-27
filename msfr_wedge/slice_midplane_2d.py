#!/usr/bin/env python3
"""Slice the real 3D wedge mesh at its angular midplane (theta=11.25 deg,
where the real inlet/outlet channel sits -- confirmed earlier from the
real mesh's own Inlet/Outlet face clustering) to get a flat 2D (r,z)
cross-section for fast 2D testing.

This is the "Flat 2D planar approximation" fallback documented at the
very start of this project (flat_2d_planar_approximation_PLAN.md): a
true FLAT 2D domain, not an axisymmetric reduction -- it loses the 1/r
geometric term an actual axisymmetric formulation would carry (Fluidity
has no native axisymmetric solver, confirmed early in this project).
Useful for fast iteration on BCs/solver/turbulence-model settings before
committing to slow 3D runs, not a physically equivalent reduced model.

Method: cut the combined (tets + tagged boundary triangles) mesh with
the theta=11.25 deg plane using vtkCutter (pyvista .slice()). Cutting a
3D tet produces a 2D polygon (the interior/bulk fill); cutting a tagged
2D boundary triangle produces a 1D line segment that inherits the
parent triangle's PhysicalTag (wall=3/inlet=1/outlet=2) -- so one slice
operation gives both the filled cross-section AND its correctly tagged
boundary in one pass. symmetry_0/symmetry_1 (tag 4/5, the theta=0 and
theta=22.5 deg planes) contribute no segments here since they don't
intersect the theta=11.25 deg midplane except trivially on the axis.

Output coordinates are reprojected onto a physically meaningful local
(r, y) basis: r = distance from the reactor's y-axis (always >=0 on this
half-plane), y = the true axial coordinate, unchanged.
"""
import numpy as np
import meshio
import pyvista as pv
import vtk

SRC = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_real_tet.msh"
OUT_MSH = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_midplane_2d.msh"
OUT_VTU = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_midplane_2d_geometry.vtu"

TAG_NAMES = {1: "inlet", 2: "outlet", 3: "wall", 4: "axis"}
FLUID_TAG = 5
FLUID_NAME = "fluid"
# Confirmed by direct inspection: the midplane cut meets symmetry_0's own
# mesh (tag 4 in the 3D source) along a real line, not just a point --
# every resulting edge has r=0.0 exactly, because the real core-bulk mesh
# extends all the way to the reactor's central axis (fluid flows through
# the centre, it isn't a solid post) and symmetry_0's triangulated patch
# happens to include that inner r=0 edge (symmetry_1, tag 5, gave zero
# edges here -- apparently an asymmetry in the real mesh's near-axis
# construction between the two cut faces, not something we rely on).
# Relabelled "axis": this is the natural boundary where our flat 2D slice
# meets the centerline, not a solid wall.

print("Reading", SRC)
mesh = meshio.read(SRC)
pts = mesh.points

tets = tris = None
tet_tags = tri_tags = None
for block, tagarray in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if block.type == "tetra":
        tets, tet_tags = block.data, tagarray
    elif block.type == "triangle":
        tris, tri_tags = block.data, tagarray

n_tet, n_tri = len(tets), len(tris)
cell_types_in = np.concatenate([
    np.full(n_tri, vtk.VTK_TRIANGLE, dtype=np.uint8),
    np.full(n_tet, vtk.VTK_TETRA, dtype=np.uint8),
])


def cellarray(conn, npts):
    n = len(conn)
    arr = np.empty((n, npts + 1), dtype=np.int64)
    arr[:, 0] = npts
    arr[:, 1:] = conn
    return arr.ravel()


ca = np.concatenate([cellarray(tris, 3), cellarray(tets, 4)])
phys = np.concatenate([tri_tags.astype(np.int32), tet_tags.astype(np.int32)])

grid = pv.UnstructuredGrid(ca, cell_types_in, pts)
grid.cell_data["PhysicalTag"] = phys

# --- cut at the real channel's angular midplane ----------------------------
theta_mid = np.radians(11.25)
n_cut = (-np.sin(theta_mid), 0.0, np.cos(theta_mid))
print(f"Cutting at theta={np.degrees(theta_mid)} deg, normal={n_cut}")

cutter = vtk.vtkCutter()
cutter.SetInputData(grid)
plane = vtk.vtkPlane()
plane.SetOrigin(0.0, 0.0, 0.0)
plane.SetNormal(*n_cut)
cutter.SetCutFunction(plane)
cutter.Update()
cut = pv.wrap(cutter.GetOutput()).cast_to_unstructured_grid()
print(f"Cut result: {cut.n_points} points, {cut.n_cells} cells")
print("Cell types present:", np.unique(cut.celltypes, return_counts=True))

cut_tags = cut.cell_data["PhysicalTag"]

# --- separate 2D fill (from tets) vs 1D boundary (from tagged triangles) ---
pts2 = cut.points
r = np.hypot(pts2[:, 0], pts2[:, 2])
y = pts2[:, 1]
local_xy = np.column_stack([r, y])  # local 2D coords: (r, axial)

fill_tris = []   # list of (i,j,k) local node index tuples
# symmetry_0/symmetry_1 (tags 4,5) shouldn't contribute real segments here
# (this plane only meets them trivially, on the axis) -- collected anyway
# and asserted empty after degenerate filtering, rather than assumed.
bnd_edges = {1: [], 2: [], 3: [], 4: [], 5: []}  # tag -> list of (i,j)

ctypes = cut.celltypes
cells_flat = cut.cells
i = 0
ci = 0
while i < len(cells_flat):
    npts_cell = cells_flat[i]
    conn = cells_flat[i + 1:i + 1 + npts_cell]
    t = ctypes[ci]
    tag = int(cut_tags[ci])
    if t == vtk.VTK_LINE:
        bnd_edges[tag].append(tuple(conn))
    elif t in (vtk.VTK_TRIANGLE, vtk.VTK_POLYGON, vtk.VTK_QUAD):
        # fan-triangulate any polygon (tet-plane cuts can be quads/pentagons)
        for k in range(1, len(conn) - 1):
            fill_tris.append((conn[0], conn[k], conn[k + 1]))
    else:
        raise RuntimeError(f"Unexpected cut cell type {t}")
    i += 1 + npts_cell
    ci += 1

print(f"Fill triangles: {len(fill_tris)}")
for tag in (1, 2, 3, 4, 5):
    print(f"  tag {tag} edges: {len(bnd_edges[tag])}")

# drop degenerate (zero-length/zero-area) cells from floating point noise
def seg_len(e):
    return np.linalg.norm(local_xy[e[0]] - local_xy[e[1]])


def tri_area2(t):
    p = local_xy[list(t)]
    return 0.5 * abs((p[1, 0] - p[0, 0]) * (p[2, 1] - p[0, 1]) - (p[2, 0] - p[0, 0]) * (p[1, 1] - p[0, 1]))


EPS = 1e-9
fill_tris = [t for t in fill_tris if tri_area2(t) > EPS]
for tag in bnd_edges:
    bnd_edges[tag] = [e for e in bnd_edges[tag] if seg_len(e) > EPS]

print("After degenerate filtering:")
print(f"  fill triangles: {len(fill_tris)}")
for tag in (1, 2, 3, 4, 5):
    print(f"  tag {tag} edges: {len(bnd_edges[tag])}")
assert not bnd_edges[5], (
    f"tag 5 (symmetry_1) produced {len(bnd_edges[5])} non-degenerate edges on the "
    f"midplane cut -- expected none (confirmed empty by direct inspection earlier)"
)

# --- prune unused points, renumber, write -----------------------------------
used = np.zeros(len(pts2), dtype=bool)
for t in fill_tris:
    for n in t:
        used[n] = True
for tag in bnd_edges:
    for e in bnd_edges[tag]:
        for n in e:
            used[n] = True

old_to_new = -np.ones(len(pts2), dtype=np.int64)
old_to_new[used] = np.arange(used.sum())
new_xy = local_xy[used]
new_pts3 = np.column_stack([new_xy, np.zeros(len(new_xy))])  # z=0 for a 2D mesh

cells = []
cell_data_phys = []
for tag in (1, 2, 3, 4):
    if not bnd_edges[tag]:
        continue
    conn = np.array([[old_to_new[n] for n in e] for e in bnd_edges[tag]])
    cells.append(meshio.CellBlock("line", conn))
    cell_data_phys.append(np.full(len(conn), tag, dtype=np.int32))

fill_conn = np.array([[old_to_new[n] for n in t] for t in fill_tris])
cells.append(meshio.CellBlock("triangle", fill_conn))
cell_data_phys.append(np.full(len(fill_conn), FLUID_TAG, dtype=np.int32))

field_data = {TAG_NAMES[t]: np.array([t, 1]) for t in (1, 2, 3, 4)}
field_data[FLUID_NAME] = np.array([FLUID_TAG, 2])

out_mesh = meshio.Mesh(
    points=new_pts3,
    cells=cells,
    cell_data={"gmsh:physical": cell_data_phys,
               "gmsh:geometrical": [np.ones(len(c.data), dtype=np.int32) for c in cells]},
    field_data=field_data,
)
out_mesh.write(OUT_MSH, file_format="gmsh22", binary=False)
print("Wrote", OUT_MSH)

vtu_cells = []
vtu_bid = []
for tag in (1, 2, 3, 4):
    if not bnd_edges[tag]:
        continue
    conn = np.array([[old_to_new[n] for n in e] for e in bnd_edges[tag]])
    vtu_cells.append(meshio.CellBlock("line", conn))
    vtu_bid.append(np.full(len(conn), tag, dtype=np.int32))
vtu_cells.append(meshio.CellBlock("triangle", fill_conn))
vtu_bid.append(np.zeros(len(fill_conn), dtype=np.int32))
vtu_mesh = meshio.Mesh(points=new_pts3, cells=vtu_cells, cell_data={"BoundaryID": vtu_bid})
vtu_mesh.write(OUT_VTU)
print("Wrote", OUT_VTU)

print(f"\nFinal 2D mesh: {len(new_pts3)} nodes, {len(fill_conn)} triangles")
for tag in (1, 2, 3, 4):
    print(f"  {TAG_NAMES[tag]} (tag {tag}): {len(bnd_edges[tag])} edges")
