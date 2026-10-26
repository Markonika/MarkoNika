#!/usr/bin/env python3
"""Extract a single 1/16 (22.5 deg) azimuthal wedge from the real full-core
MSFR mesh (msfr_fullcore_real.msh, built from idaholab/virtual_test_bed's
msfr_v2.geo), instead of re-deriving the wedge from the original .geo's
fragile hard-coded-index scripting.

History of attempts (kept here as a record, not aspirational docs):
  v1: selected hexahedra by centroid azimuth in [0,22.5) deg. Produced a
      jagged, non-planar "symmetry_0" face near the reactor's central
      axis (visible as stepping) -- azimuth is numerically meaningless
      right at r~0 (x=z=0).
  v2: tried a connectivity flood fill, blocked only by faces exactly on
      the two cut planes. This PROVED something important about the real
      mesh: there is no conformal (element-face-aligned) flat azimuthal
      cut anywhere near the axis, at ANY angle (tested 22.5, 45, 90, 180
      deg -- all leaked straight through to the full 616k-element mesh).
      INL's mesh uses an O-grid hub block near the pole whose internal
      mesh lines are not aligned to any azimuthal angle (standard
      practice to avoid degenerate elements at a pole) -- a flat cut
      through the axis necessarily slices through the interior of hub
      elements, not along existing faces.
  v3 (this version): geometrically CLIP the real mesh at the two planes
      (vtkClipDataSet via pyvista), which cuts through elements where
      needed and creates new sub-cells (tets/wedges) right at the hub,
      giving a truly flat, smooth boundary everywhere -- at the cost of
      the cut cells near the hub no longer being the exact original
      hexahedra (mesh quality there is whatever vtkClipDataSet produces,
      not independently reviewed/optimized).

Verified from the real mesh that the revolution axis is global Y, and
that each 22.5-degree bin [i*22.5, (i+1)*22.5) contains exactly ONE inlet
channel and ONE outlet channel, centered at the bin's angular midpoint
(11.25, 33.75, ...), with the two bin-boundary planes (0 and 22.5 deg)
passing through wall-only territory -- so [0, 22.5] deg is used here.
"""
import numpy as np
import meshio
import pyvista as pv
import vtk

SRC = "/home/user/MarkoNika/msfr_wedge/msfr_fullcore_real.msh"
OUT_MSH = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_real.msh"
OUT_VTU = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_real_geometry.vtu"

THETA_HI = 22.5
PLANE_TOL = 1e-6

TAG_NAMES = {1: "inlet", 2: "outlet", 3: "wall", 4: "symmetry_0", 5: "symmetry_1"}
VOL_TAG = 6
FLUID_NAME = "fluid"

HEX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
WEDGE_FACES = [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
TET_FACES = [(0, 1, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3)]
PYRAMID_FACES = [(0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)]
FACES_BY_VTKTYPE = {
    vtk.VTK_HEXAHEDRON: HEX_FACES,
    vtk.VTK_WEDGE: WEDGE_FACES,
    vtk.VTK_TETRA: TET_FACES,
    vtk.VTK_PYRAMID: PYRAMID_FACES,
}
CELLTYPE_NAME = {
    vtk.VTK_HEXAHEDRON: "hexahedron",
    vtk.VTK_WEDGE: "wedge",
    vtk.VTK_TETRA: "tetra",
    vtk.VTK_PYRAMID: "pyramid",
    vtk.VTK_TRIANGLE: "triangle",
    vtk.VTK_QUAD: "quad",
}

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

n_quad, n_hex = len(quads), len(hexes)
cell_types_in = np.concatenate([
    np.full(n_quad, vtk.VTK_QUAD, dtype=np.uint8),
    np.full(n_hex, vtk.VTK_HEXAHEDRON, dtype=np.uint8),
])


def cellarray(conn, npts):
    n = len(conn)
    arr = np.empty((n, npts + 1), dtype=np.int64)
    arr[:, 0] = npts
    arr[:, 1:] = conn
    return arr.ravel()


ca = np.concatenate([cellarray(quads, 4), cellarray(hexes, 8)])
phys_in = np.concatenate([quad_tags.astype(np.int32), hex_tags.astype(np.int32)])

grid = pv.UnstructuredGrid(ca, cell_types_in, pts)
grid.cell_data["PhysicalTag"] = phys_in

# --- geometric clip at the two wedge-bounding planes ----------------------
th = np.radians(THETA_HI)
n1 = np.array([-np.sin(th), 0.0, np.cos(th)])

print("Clipping at theta=0 plane...")
c1 = grid.clip(normal=(0, 0, 1), origin=(0, 0, 0), invert=False)  # keep z>=0
print("Clipping at theta=22.5 plane...")
c2 = c1.clip(normal=tuple(n1), origin=(0, 0, 0), invert=True)  # keep theta<=22.5

out_pts = c2.points
ctypes = c2.celltypes
tags = c2.cell_data["PhysicalTag"]
print(f"After clipping: {c2.n_cells} cells, {len(out_pts)} points")
for t in np.unique(ctypes):
    sel = ctypes == t
    print(f"  celltype {CELLTYPE_NAME.get(t, t)}: n={sel.sum()}  tags present={np.unique(tags[sel])}")

# --- separate 2D (boundary, already tagged) vs 3D (volume) cells ----------
offset = 0
cells_flat = c2.cells
tagged_2d = {}   # frozenset(node ids) -> tag
vol_cells = []   # list of (vtk_type, node_id_tuple)
i = 0
ci = 0
while i < len(cells_flat):
    n = cells_flat[i]
    conn = tuple(cells_flat[i + 1:i + 1 + n])
    t = ctypes[ci]
    if t in (vtk.VTK_TRIANGLE, vtk.VTK_QUAD):
        tagged_2d[frozenset(conn)] = int(tags[ci])
    elif t in FACES_BY_VTKTYPE:
        vol_cells.append((t, conn))
    else:
        raise RuntimeError(f"Unexpected cell type {t} after clipping")
    i += 1 + n
    ci += 1
print(f"2D tagged cells: {len(tagged_2d)}, 3D volume cells: {len(vol_cells)}")

# --- find boundary faces of the 3D cell set (face referenced exactly once) -
face_count = {}
face_example = {}
for t, conn in vol_cells:
    for f in FACES_BY_VTKTYPE[t]:
        nodes = tuple(conn[i] for i in f)
        key = frozenset(nodes)
        face_count[key] = face_count.get(key, 0) + 1
        if key not in face_example:
            face_example[key] = nodes

boundary_faces = [(k, face_example[k]) for k, c in face_count.items() if c == 1]
print(f"Boundary faces of volume cell set: {len(boundary_faces)}")

out_faces = {1: [], 2: [], 3: [], 4: [], 5: []}
n_reused = n_new0 = n_new1 = 0
for key, nodes in boundary_faces:
    if key in tagged_2d:
        out_faces[tagged_2d[key]].append(nodes)
        n_reused += 1
        continue
    p = out_pts[list(nodes)]
    d0 = np.abs(p[:, 2]).max()
    d1 = np.abs(p @ n1).max()
    if d0 < PLANE_TOL:
        out_faces[4].append(nodes)
        n_new0 += 1
    elif d1 < PLANE_TOL:
        out_faces[5].append(nodes)
        n_new1 += 1
    else:
        raise RuntimeError(f"Unclassifiable boundary face: max|z|={d0:.3e} max|plane1|={d1:.3e} pts={p.tolist()}")

print(f"Reused original tags: {n_reused}, new symmetry_0: {n_new0}, new symmetry_1: {n_new1}")
for t in (1, 2, 3):
    print(f"  {TAG_NAMES[t]}: {len(out_faces[t])} faces")

for t, name, normal_check in [(4, "symmetry_0", lambda p: np.abs(p[:, 2]).max()),
                               (5, "symmetry_1", lambda p: np.abs(p @ n1).max())]:
    worst = max((normal_check(out_pts[list(f)]) for f in out_faces[t]), default=0.0)
    print(f"Planarity check {name}: max deviation {worst:.3e} m over {len(out_faces[t])} faces")
    assert worst < PLANE_TOL

# --- prune unused points, renumber, write ----------------------------------
used = np.zeros(len(out_pts), dtype=bool)
for t, conn in vol_cells:
    for n in conn:
        used[n] = True
for t in out_faces:
    for nodes in out_faces[t]:
        for n in nodes:
            used[n] = True

old_to_new = -np.ones(len(out_pts), dtype=np.int64)
old_to_new[used] = np.arange(used.sum())
new_pts = out_pts[used]

# group volume cells by type
vol_by_type = {}
for t, conn in vol_cells:
    vol_by_type.setdefault(t, []).append([old_to_new[n] for n in conn])

def by_vertex_count(face_list):
    """Split a list of node-id tuples into tri/quad groups (vtkClipDataSet
    can leave a tagged surface -- e.g. wall -- as a mix of both)."""
    groups = {}
    for nodes in face_list:
        groups.setdefault(len(nodes), []).append(nodes)
    return groups


cells = []
cell_data_phys = []
for tag in sorted(out_faces):
    if not out_faces[tag]:
        continue
    for n_verts, group in by_vertex_count(out_faces[tag]).items():
        conn = np.array([[old_to_new[n] for n in nodes] for nodes in group])
        ctype_name = "triangle" if n_verts == 3 else "quad"
        cells.append(meshio.CellBlock(ctype_name, conn))
        cell_data_phys.append(np.full(len(conn), tag, dtype=np.int32))

for vtype, conns in vol_by_type.items():
    cells.append(meshio.CellBlock(CELLTYPE_NAME[vtype], np.array(conns)))
    cell_data_phys.append(np.full(len(conns), VOL_TAG, dtype=np.int32))

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

vtu_cells = []
vtu_bid = []
for tag in sorted(out_faces):
    if not out_faces[tag]:
        continue
    for n_verts, group in by_vertex_count(out_faces[tag]).items():
        conn = np.array([[old_to_new[n] for n in nodes] for nodes in group])
        ctype_name = "triangle" if n_verts == 3 else "quad"
        vtu_cells.append(meshio.CellBlock(ctype_name, conn))
        vtu_bid.append(np.full(len(conn), tag, dtype=np.int32))
for vtype, conns in vol_by_type.items():
    vtu_cells.append(meshio.CellBlock(CELLTYPE_NAME[vtype], np.array(conns)))
    vtu_bid.append(np.zeros(len(conns), dtype=np.int32))
vtu_mesh = meshio.Mesh(points=new_pts, cells=vtu_cells, cell_data={"BoundaryID": vtu_bid})
vtu_mesh.write(OUT_VTU)
print("Wrote", OUT_VTU)

n_vol = sum(len(c) for c in vol_by_type.values())
print(f"\nFinal wedge mesh: {len(new_pts)} nodes, {n_vol} volume cells "
      f"({', '.join(f'{CELLTYPE_NAME[t]}={len(c)}' for t, c in vol_by_type.items())})")
for t in (1, 2, 3, 4, 5):
    print(f"  {TAG_NAMES[t]} (tag {t}): {len(out_faces[t])} faces")
