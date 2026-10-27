import pyvista as pv
import numpy as np

pv.OFF_SCREEN = True

initial = pv.read("/tmp/msfr_2d_run1_archive/msfr_wedge_midplane_2d_0.vtu")
adapted = pv.read("/tmp/msfr_2d_run1_archive/msfr_wedge_midplane_2d_1.vtu")
print(f"initial: {initial.n_points} nodes, {initial.n_cells} cells")
print(f"adapted (t~1.0): {adapted.n_points} nodes, {adapted.n_cells} cells")

# full-domain side-by-side style: two separate images, same camera/extent
for mesh, name, title in [(initial, "initial", "Initial mesh (before adaptivity)"),
                           (adapted, "adapted_t1", "Adapted mesh (t~1.0s, after 9 adapts)")]:
    p = pv.Plotter(off_screen=True, window_size=(1400, 1100))
    p.add_mesh(mesh, show_edges=True, edge_color="black", line_width=0.3, color="lightsteelblue")
    p.camera_position = "xy"
    p.add_axes()
    p.add_text(f"{title}\n{mesh.n_points} nodes, {mesh.n_cells} cells", font_size=14, color="black")
    p.screenshot(f"/home/user/MarkoNika/msfr_wedge/mesh_compare_{name}.png")

# zoomed-in crops on three regions of interest: near inlet/outlet mouth,
# near the U-bend, and near the axis (bulk/calm region) -- same physical
# window for both meshes so the refinement difference is directly visible
regions = {
    "inlet_mouth": ((1.3, 1.9), (1.0, 1.4)),
    "bend": ((1.6, 2.2), (-0.3, 0.3)),
    "axis_bulk": ((0.0, 0.5), (-0.3, 0.3)),
}

for region_name, (xr, yr) in regions.items():
    for mesh, name in [(initial, "initial"), (adapted, "adapted_t1")]:
        clipped = mesh.clip_box((xr[0], xr[1], yr[0], yr[1], -1, 1), invert=False)
        p = pv.Plotter(off_screen=True, window_size=(900, 900))
        p.add_mesh(clipped, show_edges=True, edge_color="black", line_width=0.6, color="lightsteelblue")
        p.camera_position = "xy"
        p.camera.tight(padding=0.05)
        p.add_text(f"{region_name}: {name}, {clipped.n_cells} cells in view", font_size=12, color="black")
        p.screenshot(f"/home/user/MarkoNika/msfr_wedge/mesh_zoom_{region_name}_{name}.png")

print("done")
