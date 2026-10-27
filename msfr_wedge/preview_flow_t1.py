import pyvista as pv
pv.OFF_SCREEN = True
m = pv.read("/home/user/MarkoNika/msfr_wedge/msfr_wedge_midplane_2d_1.vtu")
m["VelocityMagnitude"] = m.point_data["Velocity"].copy()
import numpy as np
m["VelocityMagnitude"] = np.linalg.norm(m.point_data["Velocity"], axis=1)

p = pv.Plotter(off_screen=True, window_size=(1400, 1100))
p.add_mesh(m, scalars="VelocityMagnitude", cmap="turbo", show_edges=False, scalar_bar_args={"title": "Velocity magnitude (m/s)"})
p.camera_position = "xy"
p.add_axes()
p.screenshot("/home/user/MarkoNika/msfr_wedge/midplane_2d_velocity_t1.png")

p2 = pv.Plotter(off_screen=True, window_size=(1400, 1100))
p2.add_mesh(m, scalars="Pressure", cmap="coolwarm", show_edges=False, scalar_bar_args={"title": "Pressure (Pa)"})
p2.camera_position = "xy"
p2.add_axes()
p2.screenshot("/home/user/MarkoNika/msfr_wedge/midplane_2d_pressure_t1.png")

p3 = pv.Plotter(off_screen=True, window_size=(1400, 1100))
p3.add_mesh(m, show_edges=True, edge_color="gray", line_width=0.2, color="lightsteelblue")
p3.camera_position = "xy"
p3.add_axes()
p3.screenshot("/home/user/MarkoNika/msfr_wedge/midplane_2d_adapted_mesh_t1.png")
print("done", m.n_points, m.n_cells)
