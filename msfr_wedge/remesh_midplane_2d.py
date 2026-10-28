#!/usr/bin/env python3
"""SUPERSEDED by remesh_midplane_2d_triangle.py - kept for reference only.

This script's gmsh-based fill preserved the traced boundary loop's points
as hard geometric constraints verbatim, including some that were
extremely close together (~0.0001-0.0009 m) but not exactly coincident -
near-duplicate points vtkCutter leaves behind in slice_midplane_2d.py when
cutting through two adjacent, slightly irregular 3D boundary triangles
near the cut plane. gmsh faithfully preserved those tiny gaps as real
features and seeded sliver triangles next to them (found by actually
running Fluidity on the result: edges down to ~0.0001 m vs. the intended
0.02 m minimum, causing the CFL-limited adaptive timestep to collapse
~160x). remesh_midplane_2d_triangle.py fixes this by merging those
near-duplicate points before meshing, and uses Triangle (Shewchuk) instead
of gmsh for the actual fill, for its provable minimum-angle quality
guarantee.

Original docstring, for context:

Regenerate the 2D midplane cross-section mesh properly.

slice_midplane_2d.py produced a 2D mesh by directly reusing whatever
triangles a plane-cutter happened to leave behind when it sliced through
the 3D tet mesh -- those triangles were never independently generated or
optimized as a 2D mesh, just inherited verbatim from local 3D element
shapes. The result looked structured/fine near the axis (inherited from
the original mesh's boundary-layer resolution there) and visibly
irregular everywhere else.

This script instead extracts the REAL boundary curve (inlet/outlet/wall/
axis edges, unchanged) from that same slice, traces it into an ordered
loop, and hands it to gmsh's own 2D mesher (Frontal-Delaunay) with
smoothing/optimization passes -- a purpose-built 2D mesh of the exact
same real geometry, not a leftover slice of a 3D volume mesh.
"""
import numpy as np
import meshio
import gmsh
from collections import defaultdict

SRC = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_midplane_2d.msh"
OUT_MSH = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_midplane_2d_remeshed.msh"
OUT_VTU = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_midplane_2d_remeshed_geometry.vtu"

TAG_NAMES = {1: "inlet", 2: "outlet", 3: "wall", 4: "axis"}
FLUID_TAG = 5
FLUID_NAME = "fluid"
SNAP_TOL = 1e-7  # merge near-coincident boundary points the plane cutter
                 # left as separate instances of the same point

print("Reading", SRC)
mesh = meshio.read(SRC)
pts = mesh.points

edges = []  # (n0, n1, tag), node ids as in the source file
for block, tagarray in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if block.type == "line":
        for conn, tag in zip(block.data, tagarray):
            edges.append((int(conn[0]), int(conn[1]), int(tag)))
print(f"{len(edges)} boundary edges")

# --- snap near-coincident points, build adjacency, trace the loop ---------
used_ids = sorted(set(n for e in edges for n in e[:2]))
coords = pts[used_ids]
key = np.round(coords, 7)
_, inv = np.unique(key, axis=0, return_inverse=True)
remap = {orig: int(inv[i]) for i, orig in enumerate(used_ids)}
snapped_xy = {}
for orig, r in zip(used_ids, inv):
    snapped_xy[int(r)] = tuple(pts[orig][:2])

adj = defaultdict(list)
for a, b, t in edges:
    ra, rb = remap[a], remap[b]
    if ra == rb:
        continue
    adj[ra].append((rb, t))
    adj[rb].append((ra, t))

degrees = {n: len(v) for n, v in adj.items()}
assert all(d == 2 for d in degrees.values()), (
    f"boundary is not a clean simple loop after snapping: "
    f"degree counts {sorted(set(degrees.values()))}"
)
print(f"Clean closed loop confirmed: {len(adj)} nodes, all degree 2")

start = next(iter(adj))
loop = [start]
tags_along = []
prev = None
cur = start
while True:
    options = [(n, t) for n, t in adj[cur] if n != prev]
    if not options:
        # back at start with only the original direction available
        options = adj[cur]
    nxt, t = options[0]
    tags_along.append(t)
    if nxt == start:
        break
    loop.append(nxt)
    prev, cur = cur, nxt
assert len(loop) == len(adj), f"loop trace got {len(loop)} nodes, expected {len(adj)}"
print(f"Traced loop: {len(loop)} nodes")

# --- split into contiguous same-tag runs -----------------------------------
runs = []  # list of (tag, [node indices into `loop`, inclusive, wrapping])
run_start = 0
for i in range(1, len(tags_along)):
    if tags_along[i] != tags_along[i - 1]:
        runs.append((tags_along[i - 1], list(range(run_start, i + 1))))
        run_start = i
runs.append((tags_along[-1], list(range(run_start, len(loop))) + [0]))
# merge a possible wrap-around split (last run and first run same tag)
if len(runs) > 1 and runs[0][0] == runs[-1][0] and runs[-1][1][-1] == 0:
    first_tag, first_idx = runs[0]
    last_tag, last_idx = runs.pop()
    runs[0] = (first_tag, last_idx[:-1] + first_idx)

print(f"{len(runs)} contiguous boundary runs:")
for tag, idx in runs:
    print(f"  tag {tag} ({TAG_NAMES[tag]}): {len(idx) - 1} segments")

# --- build the gmsh model: points + lines in the exact traced order --------
gmsh.initialize()
gmsh.model.add("midplane_2d_remesh")
geo = gmsh.model.geo

point_tag = {}
for local_id in loop:
    x, y = snapped_xy[local_id]
    point_tag[local_id] = geo.addPoint(x, y, 0.0)

curve_tags_by_tag = defaultdict(list)
all_loop_curves = []
for tag, idx in runs:
    run_point_ids = [loop[i] for i in idx]
    seg_tags = []
    for a, b in zip(run_point_ids[:-1], run_point_ids[1:]):
        lt = geo.addLine(point_tag[a], point_tag[b])
        seg_tags.append(lt)
    curve_tags_by_tag[tag].extend(seg_tags)
    all_loop_curves.extend(seg_tags)

cl = geo.addCurveLoop(all_loop_curves)
surf = geo.addPlaneSurface([cl])
geo.synchronize()

for tag in (1, 2, 3, 4):
    gmsh.model.addPhysicalGroup(1, curve_tags_by_tag[tag], tag, name=TAG_NAMES[tag])
gmsh.model.addPhysicalGroup(2, [surf], FLUID_TAG, name=FLUID_NAME)

# --- mesh with a quality 2D algorithm + optimization passes ----------------
gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay
gmsh.option.setNumber("Mesh.Smoothing", 20)  # Laplacian smoothing iterations
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)

print("Meshing 2D with Frontal-Delaunay + smoothing...")
gmsh.model.mesh.generate(2)
gmsh.model.mesh.optimize("Netgen")

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(OUT_MSH)
n_nodes = len(gmsh.model.mesh.getNodes()[0])
n_tri = sum(len(e) for e in gmsh.model.mesh.getElements(2)[1])
print(f"Remeshed 2D mesh: {n_nodes} nodes, {n_tri} triangles")
gmsh.finalize()

# --- export BoundaryID-tagged VTU for preview -------------------------------
out = meshio.read(OUT_MSH)
out_cells = []
out_bid = []
for block, tagarray in zip(out.cells, out.cell_data["gmsh:physical"]):
    out_cells.append(block)
    if block.type == "triangle":
        out_bid.append(np.zeros(len(block.data), dtype=np.int32))
    else:
        out_bid.append(tagarray.astype(np.int32))
out_mesh = meshio.Mesh(points=out.points, cells=out_cells, cell_data={"BoundaryID": out_bid})
out_mesh.write(OUT_VTU)
print("Wrote", OUT_VTU)
