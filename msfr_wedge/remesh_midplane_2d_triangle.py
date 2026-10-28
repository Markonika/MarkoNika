#!/usr/bin/env python3
"""Regenerate the 2D midplane cross-section mesh with Triangle (Shewchuk),
replacing the previous gmsh-based remesh.

Root cause of the previous mesh's sliver elements (edges down to ~0.0001 m,
vs. the 0.02 m MinimumEdgeLengths floor): slice_midplane_2d.py's vtkCutter
produces independently-computed intersection points per cut 3D boundary
triangle edge. Where the original 3D wall surface mesh has two adjacent,
slightly irregular triangles near the cut plane, this leaves two points
extremely close together (~0.0001-0.0009 m apart) but not bit-identical,
so the existing 1e-7 rounding-based duplicate merge doesn't catch them.
The old remesh_midplane_2d.py then built every loop point as a fixed gmsh
Point connected by hard-constraint Lines in exact order, so gmsh faithfully
preserved these tiny inherited gaps as real geometric features and seeded
slivers next to them - this is the same bug class as commit 3882344's
"sliver tets" fix in the original 3D mesh, just reached via a different
upstream step (vtkCutter here instead of dropped zero-area triangles).

This script: traces the same boundary loop as before, but first MERGES
consecutive loop points closer than MERGE_TOL (comfortably below the real
intended ~0.025 m boundary-layer spacing, comfortably above the ~0.0001-
0.0009 m artifacts), then hands the cleaned boundary (as a PSLG, with
segment markers carrying the inlet/outlet/wall/axis tags) to Triangle's
quality-guaranteed constrained Delaunay refinement (Ruppert's algorithm:
minimum angle + maximum area bounds), instead of gmsh's Frontal-Delaunay.
"""
import numpy as np
import meshio
import triangle as tr
from collections import defaultdict

SRC = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_midplane_2d.msh"
OUT_MSH = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_midplane_2d_remeshed.msh"
OUT_VTU = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_midplane_2d_remeshed_geometry.vtu"

TAG_NAMES = {1: "inlet", 2: "outlet", 3: "wall", 4: "axis"}
FLUID_TAG = 5
FLUID_NAME = "fluid"
SNAP_TOL = 1e-7  # merge exactly-coincident-after-rounding points (as before)
MERGE_TOL = 0.002  # merge near-duplicate (not exactly coincident) consecutive
                    # loop points closer than this - the actual fix

print("Reading", SRC)
mesh = meshio.read(SRC)
pts = mesh.points

edges = []
for block, tagarray in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if block.type == "line":
        for conn, tag in zip(block.data, tagarray):
            edges.append((int(conn[0]), int(conn[1]), int(tag)))
print(f"{len(edges)} boundary edges")

# --- snap exactly-coincident points, build adjacency, trace the loop ------
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
        options = adj[cur]
    nxt, t = options[0]
    tags_along.append(t)
    if nxt == start:
        break
    loop.append(nxt)
    prev, cur = cur, nxt
assert len(loop) == len(adj), f"loop trace got {len(loop)} nodes, expected {len(adj)}"
print(f"Traced loop: {len(loop)} nodes")

# --- merge near-duplicate consecutive loop points (the actual fix) --------
# tags_along[i] is the tag of the segment from loop[i] to loop[i+1 or wrap].
# Merge loop[i+1] into loop[i] when that segment is shorter than MERGE_TOL,
# but only when the segments on both sides of the merge share the same tag
# (defensive: never silently merge across an inlet/outlet/wall/axis
# boundary even though none of the known tiny gaps do this).
n = len(loop)
coords_arr = np.array([snapped_xy[p] for p in loop])
changed = True
n_merged = 0
while changed:
    changed = False
    n = len(loop)
    coords_arr = np.array([snapped_xy[p] for p in loop])
    next_coords = np.roll(coords_arr, -1, axis=0)
    seg_len = np.linalg.norm(coords_arr - next_coords, axis=1)
    for i in range(n):
        if seg_len[i] >= MERGE_TOL:
            continue
        t_this = tags_along[i]
        t_prev = tags_along[i - 1]
        t_next = tags_along[(i + 1) % n]
        if not (t_this == t_prev == t_next):
            print(f"  skipping merge at tiny gap {seg_len[i]:.6f} m - "
                  f"crosses a tag boundary ({TAG_NAMES.get(t_prev)}/"
                  f"{TAG_NAMES.get(t_this)}/{TAG_NAMES.get(t_next)})")
            continue
        # merge loop[i+1] into loop[i]: drop loop[i+1] and its outgoing
        # segment tag, keep everything else
        j = (i + 1) % n
        del loop[j]
        del tags_along[j if j < len(tags_along) else 0]
        n_merged += 1
        changed = True
        break  # indices shifted - restart the scan
print(f"Merged {n_merged} near-duplicate point pairs (threshold {MERGE_TOL} m)")

coords_arr = np.array([snapped_xy[p] for p in loop])
next_coords = np.roll(coords_arr, -1, axis=0)
seg_len = np.linalg.norm(coords_arr - next_coords, axis=1)
print(f"Post-merge boundary segment length: min={seg_len.min():.6f} "
      f"median={np.median(seg_len):.6f} max={seg_len.max():.6f}")

# --- build Triangle's PSLG input: vertices + segments with tag markers ----
n = len(loop)
vertices = coords_arr
segments = np.array([[i, (i + 1) % n] for i in range(n)])
segment_markers = np.array(tags_along, dtype=np.int32)

pslg = {
    "vertices": vertices,
    "segments": segments,
    "segment_markers": segment_markers,
}

print("Meshing 2D with Triangle (Shewchuk) - quality Delaunay refinement "
      "(min angle 30 deg, max area 3e-4 m^2)...")
# p: respect input PSLG segments exactly: q30: Ruppert quality refinement,
# minimum angle 30 degrees: a0.0003: max triangle area bound
result = tr.triangulate(pslg, "pq30a0.0003")

out_vertices = result["vertices"]
out_triangles = result["triangles"]
print(f"Triangle mesh: {len(out_vertices)} nodes, {len(out_triangles)} triangles")

# recover boundary edges + their tags from the refined segments (Triangle
# subdivides long input segments during refinement but keeps their marker)
out_segments = result["segments"]
out_segment_markers = result["segment_markers"].reshape(-1)

# --- assemble a meshio Mesh and write gmsh v2.2 + tagged VTU --------------
points3d = np.hstack([out_vertices, np.zeros((len(out_vertices), 1))])

line_cells = out_segments.astype(np.int64)
line_tags = out_segment_markers.astype(np.int32)
tri_cells = out_triangles.astype(np.int64)
tri_tags = np.full(len(tri_cells), FLUID_TAG, dtype=np.int32)

cells = [
    meshio.CellBlock("line", line_cells),
    meshio.CellBlock("triangle", tri_cells),
]
cell_data = {"gmsh:physical": [line_tags, tri_tags],
             "gmsh:geometrical": [line_tags, tri_tags]}

field_data = {TAG_NAMES[t]: np.array([t, 1]) for t in TAG_NAMES}
field_data[FLUID_NAME] = np.array([FLUID_TAG, 2])

out_mesh = meshio.Mesh(points=points3d, cells=cells, cell_data=cell_data,
                        field_data=field_data)
out_mesh.write(OUT_MSH, file_format="gmsh22", binary=False)
print("Wrote", OUT_MSH)

vtu_bid = [line_tags, np.zeros(len(tri_cells), dtype=np.int32)]
vtu_mesh = meshio.Mesh(points=points3d, cells=cells,
                        cell_data={"BoundaryID": vtu_bid})
vtu_mesh.write(OUT_VTU)
print("Wrote", OUT_VTU)
