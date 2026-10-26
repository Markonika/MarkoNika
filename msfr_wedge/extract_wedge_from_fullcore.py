#!/usr/bin/env python3
"""Extract a single 1/16 (22.5 deg) azimuthal wedge from the real full-core
MSFR mesh (msfr_fullcore_real.msh, built from idaholab/virtual_test_bed's
msfr_v2.geo), instead of re-deriving the wedge from the original .geo's
fragile hard-coded-index scripting.

Why slice the finished mesh rather than edit the .geo: msfr_v2.geo builds
the full 360 degree assembly with every Curve Loop/Surface/Volume index
hand-computed relative to prior entities across all 16 wedges -- there is
no clean single-wedge entry point to call. Slicing the real, already-built
mesh instead guarantees we keep the exact real geometry (true neck -> elbow
-> extrusion channel shape, real curved wall arcs) with zero re-derivation
risk.

Verified from the real mesh (not assumed from reading the script) that the
revolution axis is global Y, and that each 22.5-degree bin [i*22.5, (i+1)*22.5)
contains exactly ONE inlet channel and ONE outlet channel, centered at the
bin's angular midpoint (11.25, 33.75, ...), with clean wall-only cuts at the
bin boundaries (0, 22.5, 45, ...). So the wedge [0, 22.5) deg is used here.

Algorithm:
  1. Select all hexahedra whose centroid azimuth (atan2(z,x) mod 360) falls
     in [0, 22.5).
  2. Build a face-adjacency count over the SELECTED hexahedra's 6 quad
     faces each (by frozenset of global node ids). Any face referenced by
     exactly one selected hexahedron is on the subset's boundary.
  3. Classify each boundary face: if it matches an existing tagged quad
     element from the original mesh (Inlet=1/Outlet=2/Wall=3), keep that
     tag. Otherwise it is a newly-exposed cut face -- tag it symmetry_0
     (near theta=0) or symmetry_1 (near theta=22.5) by its centroid angle.
  4. Prune unused points, renumber, write a gmsh v2.2 .msh (for Fluidity)
     and a BoundaryID-tagged .vtu (for ParaView).
"""
import numpy as np
import meshio

SRC = "/home/user/MarkoNika/msfr_wedge/msfr_fullcore_real.msh"
OUT_MSH = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_real.msh"
OUT_VTU = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_real_geometry.vtu"

THETA_LO, THETA_HI = 0.0, 22.5
TOL = 1e-6

HEX_FACES = [
    (0, 1, 2, 3),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
]

TAG_NAMES = {1: "inlet", 2: "outlet", 3: "wall", 4: "symmetry_0", 5: "symmetry_1"}
VOL_TAG = 6
FLUID_NAME = "fluid"

print("Reading", SRC)
mesh = meshio.read(SRC)
pts = mesh.points

hexes = quads = None
hex_tags = quad_tags = None
for block, tagarray in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if block.type == "hexahedron":
        hexes, hex_tags = block.data, tagarray
    elif block.type == "quad":
        quads, quad_tags = block.data, tagarray
assert hexes is not None and quads is not None

# --- 1. select hexahedra in the wedge ---------------------------------
hex_centroids = pts[hexes].mean(axis=1)
hex_theta = np.degrees(np.arctan2(hex_centroids[:, 2], hex_centroids[:, 0])) % 360
sel = (hex_theta >= THETA_LO - TOL) & (hex_theta < THETA_HI - TOL)
sel_hexes = hexes[sel]
print(f"Selected {len(sel_hexes)} / {len(hexes)} hexahedra in [{THETA_LO},{THETA_HI}) deg")

# --- lookup for already-tagged quads (Inlet/Outlet/Wall) --------------
orig_quad_tag = {}
for q, t in zip(quads, quad_tags):
    orig_quad_tag[frozenset(q.tolist())] = (tuple(q.tolist()), int(t))

# --- 2. face adjacency count over selected hexahedra -------------------
face_count = {}
face_example = {}
for hx in sel_hexes:
    for f in HEX_FACES:
        nodes = tuple(hx[i] for i in f)
        key = frozenset(nodes)
        face_count[key] = face_count.get(key, 0) + 1
        if key not in face_example:
            face_example[key] = nodes

boundary_faces = [face_example[k] for k, c in face_count.items() if c == 1]
print(f"Boundary faces on the wedge subset: {len(boundary_faces)}")

# --- 3. classify boundary faces ----------------------------------------
out_quads = {1: [], 2: [], 3: [], 4: [], 5: []}
n_reused = 0
n_new_sym0 = 0
n_new_sym1 = 0
n_unclassified = 0
for nodes in boundary_faces:
    key = frozenset(nodes)
    if key in orig_quad_tag:
        orig_nodes, tag = orig_quad_tag[key]
        out_quads[tag].append(orig_nodes)
        n_reused += 1
    else:
        c = pts[list(nodes)].mean(axis=0)
        theta = np.degrees(np.arctan2(c[2], c[0])) % 360
        # distance (in angle) to each candidate symmetry plane
        d0 = min(theta, 360 - theta)
        d1 = abs(theta - THETA_HI)
        if d0 <= d1:
            out_quads[4].append(nodes)
            n_new_sym0 += 1
        else:
            out_quads[5].append(nodes)
            n_new_sym1 += 1

print(f"Reused original tags (inlet/outlet/wall): {n_reused}")
print(f"New symmetry_0 (theta~0) faces: {n_new_sym0}")
print(f"New symmetry_1 (theta~22.5) faces: {n_new_sym1}")
for t in (1, 2, 3):
    print(f"  {TAG_NAMES[t]}: {len(out_quads[t])} quads")

# --- 4. prune unused points, renumber, write ----------------------------
used = np.zeros(len(pts), dtype=bool)
used[sel_hexes.ravel()] = True
for t in out_quads:
    for nodes in out_quads[t]:
        for n in nodes:
            used[n] = True

old_to_new = -np.ones(len(pts), dtype=np.int64)
old_to_new[used] = np.arange(used.sum())
new_pts = pts[used]

new_hex = old_to_new[sel_hexes]
new_quads_by_tag = {t: old_to_new[np.array(v, dtype=np.int64)] for t, v in out_quads.items() if v}

cells = []
cell_data_phys = []
for t in sorted(new_quads_by_tag):
    cells.append(meshio.CellBlock("quad", new_quads_by_tag[t]))
    cell_data_phys.append(np.full(len(new_quads_by_tag[t]), t, dtype=np.int32))
cells.append(meshio.CellBlock("hexahedron", new_hex))
cell_data_phys.append(np.full(len(new_hex), VOL_TAG, dtype=np.int32))

field_data = {TAG_NAMES[t]: np.array([t, 2]) for t in (1, 2, 3, 4, 5)}
field_data[FLUID_NAME] = np.array([VOL_TAG, 3])

out_mesh = meshio.Mesh(
    points=new_pts,
    cells=cells,
    cell_data={"gmsh:physical": cell_data_phys,
               "gmsh:geometrical": [np.ones(len(c.data), dtype=np.int32) for c in cells]},
    field_data=field_data,
)
out_mesh.write(OUT_MSH, file_format="gmsh22", binary=False)
print("Wrote", OUT_MSH)

# also a BoundaryID-tagged VTU for ParaView (surface quads + volume hexes,
# interior hexes get BoundaryID=0 so a Threshold filter can isolate the
# boundary surfaces, matching the earlier full-core export convention).
vtu_bid = []
for t in sorted(new_quads_by_tag):
    vtu_bid.append(np.full(len(new_quads_by_tag[t]), t, dtype=np.int32))
vtu_bid.append(np.zeros(len(new_hex), dtype=np.int32))
vtu_cells = [meshio.CellBlock("quad", new_quads_by_tag[t]) for t in sorted(new_quads_by_tag)]
vtu_cells.append(meshio.CellBlock("hexahedron", new_hex))
vtu_mesh = meshio.Mesh(points=new_pts, cells=vtu_cells, cell_data={"BoundaryID": vtu_bid})
vtu_mesh.write(OUT_VTU)
print("Wrote", OUT_VTU)

print(f"\nFinal wedge mesh: {len(new_pts)} nodes, {len(new_hex)} hexahedra")
for t in (1, 2, 3, 4, 5):
    n = len(new_quads_by_tag.get(t, []))
    print(f"  {TAG_NAMES[t]} (tag {t}): {n} quads")
