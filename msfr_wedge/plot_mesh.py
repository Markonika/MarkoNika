import gmsh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

gmsh.initialize()
gmsh.open("/home/user/MarkoNika/msfr_wedge/msfr_wedge.msh")

node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
coords = {t: node_coords[3 * i:3 * i + 3] for i, t in enumerate(node_tags)}

fig = plt.figure(figsize=(14, 7))
ax3d = fig.add_subplot(121, projection="3d")
ax2d = fig.add_subplot(122)

colors = {"inlet": "red", "outlet": "blue", "wall": "lightgray",
          "symmetry_0": "orange", "symmetry_1": "gold"}

for dim, tag in gmsh.model.getPhysicalGroups(2):
    name = gmsh.model.getPhysicalName(dim, tag)
    ents = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
    tris = []
    for ent in ents:
        etypes, etags, enodes = gmsh.model.mesh.getElements(2, int(ent))
        for et, nodes in zip(etypes, enodes):
            if et == 2:  # 3-node triangle
                nodes = nodes.reshape(-1, 3)
                for tri in nodes:
                    tris.append([coords[n] for n in tri])
    if tris:
        pc = Poly3DCollection(tris, facecolor=colors.get(name, "green"),
                               edgecolor="k", linewidths=0.1, alpha=0.9)
        ax3d.add_collection3d(pc)

all_pts = np.array(list(coords.values()))
ax3d.set_xlim(all_pts[:, 0].min(), all_pts[:, 0].max())
ax3d.set_ylim(all_pts[:, 1].min(), all_pts[:, 1].max())
ax3d.set_zlim(all_pts[:, 2].min(), all_pts[:, 2].max())
ax3d.set_box_aspect([1, 1, 1])
ax3d.view_init(elev=15, azim=-60)
ax3d.set_title("3D wedge mesh (red=inlet, blue=outlet, gray=wall)")

# 2D (r,z) view of the symmetry_0 face (the original profile)
tag_sym0 = [t for d, t in gmsh.model.getPhysicalGroups(2) if gmsh.model.getPhysicalName(d, t) == "symmetry_0"][0]
ents = gmsh.model.getEntitiesForPhysicalGroup(2, tag_sym0)
for ent in ents:
    etypes, etags, enodes = gmsh.model.mesh.getElements(2, int(ent))
    for et, nodes in zip(etypes, enodes):
        if et == 2:
            nodes = nodes.reshape(-1, 3)
            for tri in nodes:
                pts = np.array([coords[n] for n in tri])
                r = np.hypot(pts[:, 1], pts[:, 2])
                z = pts[:, 0]
                ax2d.fill(z, r, facecolor="lightsteelblue", edgecolor="k", linewidth=0.3)
ax2d.set_xlabel("z (axial, m)")
ax2d.set_ylabel("r (radial, m)")
ax2d.set_title("Profile cross-section (r vs z)")
ax2d.set_aspect("equal")
ax2d.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("/home/user/MarkoNika/msfr_wedge/mesh_preview.png", dpi=150)
print("saved mesh_preview.png")
gmsh.finalize()
