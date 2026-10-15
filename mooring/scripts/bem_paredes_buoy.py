#!/usr/bin/env python3
"""Linear BEM solution (Capytaine) for the hull of the Paredes (2016) buoy: added mass, radiation damping and wave excitation (surge, heave, pitch)
over a range of periods in 0.9 m of water. This is *computed* data (not a measurement); see docs/bem_paredes.md.

Geometry (from the thesis data used in docs/paredes_validation.md): vertical cylinder D = 0.515 m (C33 = rho g pi D^2/4 = 2043.5 N/m), draft 0.176 m,
centre of gravity 0.1002 m below the waterline (examples/paredes: cg_ref 0.7998 m with the free surface at 0.9 m), depth 0.9 m, rho = 1000.
Output: JSON table with, per period, A and B (6x6, only the surge-heave-pitch block is non-zero) and the excitation w = |F|, delta = -arg(F) in the
convention F(t) = w a sin(omega t + delta) for an elevation a sin(omega t) at x = 0 (waves along +x), the convention of the mooring code.
usage: bem_paredes_buoy.py OUT.json [mesh_resolution_factor]"""
import json, sys
import numpy as np
import capytaine as cpt
from capytaine.bem.airy_waves import froude_krylov_force

out = sys.argv[1]
fac = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
rho, g, depth = 1000.0, 9.81, 0.9
D = np.sqrt(4.0 * 2043.494 / (rho * g * np.pi))
R, draft, zcg = D / 2.0, 0.176, -0.1002
# The mesh is a taller cylinder clipped at the waterline: a panel lying in the free-surface plane (the top cap) would spoil the Green-function solution.
top = 0.12
res = (max(2, int(6 * fac)), max(8, int(48 * fac)), max(4, int(round(14 * fac * (draft + top) / draft))))
mesh = cpt.mesh_vertical_cylinder(length=draft + top, radius=R, center=(0, 0, -draft + (draft + top) / 2.0), resolution=res)
body = cpt.FloatingBody(mesh, dofs=cpt.rigid_body_dofs(rotation_center=(0.0, 0.0, zcg))).immersed_part(free_surface=0.0)
print("D = %.4f m, mesh %s -> %d wetted panels, displaced volume %.5f m^3 (cylinder %.5f)" % (D, res, body.mesh.nb_faces, body.volume, np.pi * R ** 2 * draft), file=sys.stderr)

periods = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.3, 2.6, 3.0, 3.5, 4.0, 5.0]
solver = cpt.BEMSolver()
names = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]
rows = []
for T in periods:
    w = 2.0 * np.pi / T
    A = np.zeros((6, 6)); B = np.zeros((6, 6)); F = np.zeros(6, dtype=complex)
    for i, di in enumerate(names):
        res_i = solver.solve(cpt.RadiationProblem(body=body, radiating_dof=di, omega=w, water_depth=depth, rho=rho, g=g))
        for j, dj in enumerate(names):
            A[j, i] = res_i.added_masses[dj]; B[j, i] = res_i.radiation_dampings[dj]
    dpb = cpt.DiffractionProblem(body=body, omega=w, wave_direction=0.0, water_depth=depth, rho=rho, g=g)
    dres = solver.solve(dpb)
    fk = froude_krylov_force(dpb)                    # DiffractionResult.forces holds only the diffraction part: the excitation is FK + diffraction
    for j, dj in enumerate(names):
        F[j] = dres.forces[dj] + fk[dj]
    # Capytaine: e^{-i w t}; an elevation a sin(w t - k x) corresponds to the complex amplitude i a  ->  F(t) = a |F| sin(w t - arg F)
    rows.append({"period_s": T, "A": A.tolist(), "B": B.tolist(), "w": np.abs(F).tolist(), "delta": (-np.angle(F)).tolist()})
    print("T %.2f  A11 %.3f A33 %.3f A55 %.4f A15 %.4f | B11 %.3f B33 %.3f B55 %.5f B15 %.4f | w1 %.1f d1 %.3f w3 %.1f d3 %.3f w5 %.2f d5 %.3f" % (
        T, A[0, 0], A[2, 2], A[4, 4], A[0, 4], B[0, 0], B[2, 2], B[4, 4], B[0, 4], np.abs(F[0]), -np.angle(F[0]), np.abs(F[2]), -np.angle(F[2]), np.abs(F[4]), -np.angle(F[4])), file=sys.stderr)
json.dump({"_comment": "Capytaine BEM, Paredes hull, depth 0.9 m, CG 0.1002 m below the waterline; see scripts/bem_paredes_buoy.py", "rows": rows}, open(out, "w"), indent=1)
