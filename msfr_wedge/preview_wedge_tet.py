import pyvista as pv
pv.OFF_SCREEN = True
mesh = pv.read("/home/user/MarkoNika/msfr_wedge/msfr_wedge_real_tet_geometry.vtu")
print(mesh)

colors = {1: "red", 2: "blue", 3: "lightgray", 4: "orange", 5: "gold"}
p = pv.Plotter(off_screen=True, window_size=(1400, 1100))
for tag, name in [(3, "wall"), (1, "inlet"), (2, "outlet"), (4, "symmetry_0"), (5, "symmetry_1")]:
    part = mesh.threshold([tag, tag], scalars="BoundaryID")
    if part.n_cells == 0:
        continue
    opacity = 1.0 if tag in (1, 2) else (0.25 if tag == 3 else 0.5)
    p.add_mesh(part.extract_surface(), color=colors[tag], opacity=opacity,
               show_edges=(tag in (1, 2)), label=name)
p.add_legend()
p.camera_position = "iso"
p.add_axes()
p.screenshot("/home/user/MarkoNika/msfr_wedge/wedge_tet_preview_tagged.png")

clean = mesh.threshold([1, 5], scalars="BoundaryID")
p2 = pv.Plotter(off_screen=True, window_size=(1400, 1100))
p2.add_mesh(clean.extract_surface(), scalars="BoundaryID", cmap="tab10", show_edges=True, line_width=0.3)
p2.camera_position = "xy"
p2.add_axes()
p2.screenshot("/home/user/MarkoNika/msfr_wedge/wedge_tet_preview_top.png")

# cut the volume open to see tet quality/structure in the interior
vol_only = mesh.threshold([0, 0], scalars="BoundaryID")
clip = vol_only.clip(normal=(0, 1, 0))
p3 = pv.Plotter(off_screen=True, window_size=(1400, 1100))
p3.add_mesh(clip, show_edges=True, color="lightblue", edge_color="gray", line_width=0.3)
p3.camera_position = "iso"
p3.add_axes()
p3.screenshot("/home/user/MarkoNika/msfr_wedge/wedge_tet_preview_interior.png")
print("done")
