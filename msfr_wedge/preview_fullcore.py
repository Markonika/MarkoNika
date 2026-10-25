#!/usr/bin/env python3
import pyvista as pv

pv.OFF_SCREEN = True
mesh = pv.read("/home/user/MarkoNika/msfr_wedge/msfr_fullcore_real_geometry.vtu")
print(mesh)
surf = mesh.extract_surface()

plotter = pv.Plotter(off_screen=True, window_size=(1400, 1100))
plotter.add_mesh(surf, scalars="BoundaryID", cmap="tab10", show_edges=False)
plotter.camera_position = "iso"
plotter.add_axes()
plotter.screenshot("/home/user/MarkoNika/msfr_wedge/fullcore_preview_iso.png")

plotter2 = pv.Plotter(off_screen=True, window_size=(1400, 1100))
plotter2.add_mesh(surf, scalars="BoundaryID", cmap="tab10", show_edges=False)
plotter2.camera_position = "xy"
plotter2.add_axes()
plotter2.screenshot("/home/user/MarkoNika/msfr_wedge/fullcore_preview_top.png")
print("done")
