import pyvista as pv
pv.OFF_SCREEN = True
mesh = pv.read("/home/user/MarkoNika/msfr_wedge/msfr_wedge_midplane_2d_remeshed_geometry.vtu")
print(mesh)

fill = mesh.threshold([0, 0], scalars="BoundaryID")
p = pv.Plotter(off_screen=True, window_size=(1400, 1100))
p.add_mesh(fill, color="lightsteelblue", show_edges=True, edge_color="gray", line_width=0.3)

colors = {1: "red", 2: "blue", 3: "black", 4: "green"}
names = {1: "inlet", 2: "outlet", 3: "wall", 4: "axis"}
for tag in (1, 2, 3, 4):
    part = mesh.threshold([tag, tag], scalars="BoundaryID")
    if part.n_cells == 0:
        continue
    p.add_mesh(part, color=colors[tag], line_width=4, label=names[tag])
p.add_legend()
p.camera_position = "xy"
p.add_axes()
p.screenshot("/home/user/MarkoNika/msfr_wedge/midplane_2d_remeshed_preview.png")
print("done")
