import meshio
import numpy as np
import matplotlib.pyplot as plt

m = meshio.read("/home/user/MarkoNika/msfr_wedge/msfr_wedge_6.vtu")
pts = m.points
vel = m.point_data["Velocity"]

# project onto the symmetry_0 plane (z_gmsh ~ 0, i.e. points near theta=0)
# the wedge is thin (4 deg), so selecting points with the smallest
# |gmsh z| gives (approximately) the theta=0 face
mask = pts[:, 2] < 1e-9

x = pts[mask, 0]            # axial
r = pts[mask, 1]            # radial (theta~0, so y~r)
u_ax = vel[mask, 0]         # axial velocity
u_rad = vel[mask, 1]        # radial-ish velocity (y-component)
speed = np.linalg.norm(vel[mask], axis=1)

fig, axs = plt.subplots(1, 2, figsize=(14, 6))

sc = axs[0].scatter(x, r, c=speed, cmap="turbo", s=14)
axs[0].set_xlabel("z (axial, m)")
axs[0].set_ylabel("r (radial, m)")
axs[0].set_title(f"Velocity magnitude [m/s] at t={m.field_data.get('Time', ['?'])}\n(symmetry_0 slice, t~0.84s into transient)")
axs[0].set_aspect("equal")
plt.colorbar(sc, ax=axs[0])

step = max(1, len(x) // 400)
axs[1].quiver(x[::step], r[::step], u_ax[::step], u_rad[::step],
              speed[::step], cmap="turbo", scale=2.0, width=0.003)
axs[1].set_xlabel("z (axial, m)")
axs[1].set_ylabel("r (radial, m)")
axs[1].set_title("Velocity vectors (axial, radial)")
axs[1].set_aspect("equal")

plt.tight_layout()
plt.savefig("/home/user/MarkoNika/msfr_wedge/flow_preview.png", dpi=150)
print("saved flow_preview.png")
print("speed range:", speed.min(), speed.max())
print("n points on slice:", mask.sum())
