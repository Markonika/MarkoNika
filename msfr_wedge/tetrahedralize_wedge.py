#!/usr/bin/env python3
"""Rebuild the real MSFR 1/16 wedge as an all-tetrahedra mesh.

Keeps the exact boundary surface we already verified (real wall/inlet/
outlet geometry from idaholab/virtual_test_bed, symmetry_0/symmetry_1
flat to 1e-16 m from the geometric clip) and only replaces the interior:
quad boundary faces are split into triangles, fed into gmsh as 5 discrete
surfaces (one per physical tag), stitched into a closed shell, and the
enclosed volume is filled with gmsh's own tetrahedral mesher. This avoids
naively splitting the existing hexahedra/wedges into tets, which would
risk non-conforming faces between neighbouring cells (adjacent hexes can
pick inconsistent face diagonals unless you're careful) -- re-meshing the
interior from a fixed boundary sidesteps that entirely.
"""
import numpy as np
import meshio
import gmsh

SRC = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_real.msh"
OUT_MSH = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_real_tet.msh"
OUT_VTU = "/home/user/MarkoNika/msfr_wedge/msfr_wedge_real_tet_geometry.vtu"

TAG_NAMES = {1: "inlet", 2: "outlet", 3: "wall", 4: "symmetry_0", 5: "symmetry_1"}

print("Reading", SRC)
mesh = meshio.read(SRC)
pts = mesh.points

def tri_area(tri):
    p = pts[list(tri)]
    return 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))


AREA_EPS = 1e-12

# gather all boundary (quad/triangle) elements per physical tag, splitting
# quads into two triangles each. Some wall quads in the real source mesh
# are degenerate (two corners coincide -- a boundary-layer cell pinched
# to zero thickness at a sharp corner, 378 of them here), which produces
# a zero-area "sliver" triangle on one of the two diagonal splits; that
# sliver creates a non-manifold edge once boundary patches are stitched
# together (an edge shared by far more than 2 triangles), which gmsh's
# createGeometry() rejects outright. Drop zero-area triangles instead.
tris_by_tag = {t: [] for t in TAG_NAMES}
n_dropped = 0
for block, tagarray in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if block.type == "triangle":
        for conn, tag in zip(block.data, tagarray):
            tri = tuple(conn)
            if tri_area(tri) > AREA_EPS:
                tris_by_tag[int(tag)].append(tri)
            else:
                n_dropped += 1
    elif block.type == "quad":
        for conn, tag in zip(block.data, tagarray):
            a, b, c, d = conn
            for tri in [(a, b, c), (a, c, d)]:
                if tri_area(tri) > AREA_EPS:
                    tris_by_tag[int(tag)].append(tri)
                else:
                    n_dropped += 1

print(f"Dropped {n_dropped} degenerate (zero-area) triangles")
for t in TAG_NAMES:
    print(f"  tag {t} ({TAG_NAMES[t]}): {len(tris_by_tag[t])} triangles")

gmsh.initialize()
gmsh.model.add("wedge_tet")

node_tag_counter = 1
surf_tags = []
for tag in sorted(TAG_NAMES):
    tris = tris_by_tag[tag]
    if not tris:
        continue
    gmsh.model.addDiscreteEntity(2, tag)
    surf_tags.append(tag)

    used = sorted(set(n for tri in tris for n in tri))
    local_id = {n: node_tag_counter + i for i, n in enumerate(used)}
    node_tag_counter += len(used)

    node_tags = np.array([local_id[n] for n in used], dtype=np.int64)
    coords = pts[used].ravel()
    gmsh.model.mesh.addNodes(2, tag, node_tags, coords)

    tri_conn = np.array([[local_id[n] for n in tri] for tri in tris], dtype=np.int64)
    elem_tags = np.arange(1, len(tris) + 1, dtype=np.int64) + (tag * 10_000_000)
    gmsh.model.mesh.addElements(2, tag, [2], [elem_tags], [tri_conn.ravel()])

print("Merging coincident boundary nodes across patches...")
gmsh.model.mesh.removeDuplicateNodes()

sl = gmsh.model.geo.addSurfaceLoop(surf_tags)
vol = gmsh.model.geo.addVolume([sl])
gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(3, [vol], 6, name="fluid")
for tag in surf_tags:
    gmsh.model.addPhysicalGroup(2, [tag], tag, name=TAG_NAMES[tag])

print("Meshing volume with tetrahedra...")
gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay
gmsh.model.mesh.generate(3)

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(OUT_MSH)

n_nodes = len(gmsh.model.mesh.getNodes()[0])
n_tet = sum(len(e) for e in gmsh.model.mesh.getElements(3)[1])
print(f"Tet mesh written: {n_nodes} nodes, {n_tet} tetrahedra")
gmsh.finalize()

# --- export a BoundaryID-tagged VTU for ParaView ---------------------------
out = meshio.read(OUT_MSH)
out_cells = []
out_bid = []
for block, tagarray in zip(out.cells, out.cell_data["gmsh:physical"]):
    out_cells.append(block)
    if block.type == "tetra":
        out_bid.append(np.zeros(len(block.data), dtype=np.int32))
    else:
        out_bid.append(tagarray.astype(np.int32))
out_mesh = meshio.Mesh(points=out.points, cells=out_cells, cell_data={"BoundaryID": out_bid})
out_mesh.write(OUT_VTU)
print("Wrote", OUT_VTU)
