#!/usr/bin/env python3
import pyvista as pv

pv.OFF_SCREEN = True
mesh = pv.read("/home/user/MarkoNika/msfr_wedge/msfr_wedge_real_geometry.vtu")
print(mesh)
surf = mesh.extract_surface()

p = pv.Plotter(off_screen=True, window_size=(1400, 1100))
p.add_mesh(surf, scalars="BoundaryID", cmap="tab10", show_edges=True, edge_color="gray", line_width=0.3)
p.camera_position = "iso"
p.add_axes()
p.screenshot("/home/user/MarkoNika/msfr_wedge/wedge_real_preview_iso.png")

# only the true boundary (BoundaryID>0), no interior hex faces -- avoids
# z-fighting, matches the ParaView Threshold trick we told the user about
clean = mesh.threshold([1, 5], scalars="BoundaryID")
p2 = pv.Plotter(off_screen=True, window_size=(1400, 1100))
p2.add_mesh(clean.extract_surface(), scalars="BoundaryID", cmap="tab10", show_edges=True, edge_color="black", line_width=0.4)
p2.camera_position = "iso"
p2.add_axes()
p2.screenshot("/home/user/MarkoNika/msfr_wedge/wedge_real_preview_clean.png")

p3 = pv.Plotter(off_screen=True, window_size=(1400, 1100))
p3.add_mesh(clean.extract_surface(), scalars="BoundaryID", cmap="tab10", show_edges=True, edge_color="black", line_width=0.4)
p3.camera_position = "xy"
p3.camera.azimuth = 0
p3.add_axes()
p3.screenshot("/home/user/MarkoNika/msfr_wedge/wedge_real_preview_top.png")
print("done")
