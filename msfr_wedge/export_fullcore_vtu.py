#!/usr/bin/env python3
"""Export the real full-core MSFR mesh (from idaholab/virtual_test_bed's
msfr_v2.geo) to a VTU with a BoundaryID cell-data array, for viewing in
ParaView. Physical groups per the .geo: Inlet=1, Outlet=2, Wall=3 (surfaces),
Fluid=4 (volume).
"""
import meshio

mesh = meshio.read("/home/user/MarkoNika/msfr_wedge/msfr_fullcore_real.msh")

cells = []
boundary_id = []
for cb in mesh.cells:
    tags = mesh.cell_data["gmsh:physical"][len(cells)] if False else None

# meshio keeps cell blocks in mesh.cells and matching gmsh:physical in
# mesh.cell_data["gmsh:physical"], aligned by block index.
phys = mesh.cell_data["gmsh:physical"]
out_cells = []
out_bid = []
for block, tagarray in zip(mesh.cells, phys):
    out_cells.append(block)
    out_bid.append(tagarray)

out_mesh = meshio.Mesh(
    points=mesh.points,
    cells=out_cells,
    cell_data={"BoundaryID": out_bid},
)
out_path = "/home/user/MarkoNika/msfr_wedge/msfr_fullcore_real_geometry.vtu"
out_mesh.write(out_path)

n_hex = sum(len(b) for b in mesh.cells if b.type == "hexahedron")
n_quad = sum(len(b) for b in mesh.cells if b.type == "quad")
print(f"nodes={len(mesh.points)} hex={n_hex} quad(surface)={n_quad}")
print(f"wrote {out_path}")
