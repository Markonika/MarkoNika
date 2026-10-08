#!/usr/bin/env python3
"""Generate examples/paredes/*.json from the data of Paredes (2016), 'Study of Mooring Systems for Offshore ...'.
Every number is annotated with its source; items marked ASSUMPTION are not in the thesis (see docs/paredes_validation.md)."""
import json, math, os
here = os.path.dirname(os.path.abspath(__file__)); out = os.path.join(os.path.dirname(here), "examples", "paredes")
g, rho = 9.81, 1000.0

# --- buoy (Table 5.3, 5.4, 5.12-5.14, sec. 5.3) ---
m_b, D_b, H_b, I_cg, y_cg = 35.50, 0.515, 0.400, 0.87, 0.0758          # Table 5.3
draft_free = 0.176                                                      # Table 5.4
depth = 0.900                                                           # Annex A.1
C33 = rho * g * math.pi * D_b**2 / 4                                    # waterplane stiffness (task brief; thesis Kt = 2039.5)
C55 = 37.76                                                             # Table 5.14, free buoy (computed by the thesis)
r_fair = D_b / 2 + 0.015                                                # fairlead 0.015 m from the hull (sec. 5.3)
z_fair_body = (H_b - 0.220) - y_cg                                      # 0.220 m below the top of the buoy (sec. 5.3), relative to the CG
cg_ref = [0.0, 0.0, depth - draft_free + y_cg]
# potential-theory added mass and radiation damping per configuration (Tables 5.12-5.14: mab, brb)
hyd = {
  "free": dict(a11=18.28, a33=26.92, a55=0.283, b11=0.0196, b33=37.92, b55=0.0493),   # surge value unused (no restoring force)
  "con1": dict(a11=18.28, a33=26.96, a55=0.285, b11=0.0196, b33=37.94, b55=0.0605),
  "con2": dict(a11=20.14, a33=27.11, a55=0.311, b11=0.0182, b33=35.91, b55=0.1019),
  "cat":  dict(a11=18.73, a33=27.06, a55=0.291, b11=0.0167, b33=36.04, b55=0.0655)}
# floaters / clumpweights (Tables 5.7, 5.8)
floaters = [(0.252, 9.97), (0.276, 10.05), (0.232, 9.96)]
clumps = [(1.115, 9.968), (1.125, 10.051), (1.113, 9.959)]
D_floater = 0.135

def body(k):
    h = hyd[k]
    return {"mass_kg": m_b, "inertia_diag_kg_m2": [I_cg, I_cg, 1.18],     # yaw inertia: ASSUMPTION (solid-cylinder estimate)
            "cg_ref_m": cg_ref,
            "A_diag": [h["a11"], h["a11"], h["a33"], h["a55"], h["a55"], 0.0],
            "B_diag": [h["b11"], h["b11"], h["b33"], h["b55"], h["b55"], 0.0],
            "C_diag": [0.0, 0.0, round(C33, 3), C55, C55, 0.0], "Dq": [0, 0, 0, 0, 0, 0], "F0": [0, 0, 0, 0, 0, 0]}

def leg_geometry(k_leg, horiz):
    ang = [120.0, 0.0, -120.0][k_leg]                                   # leg 2 leeward (+x = wave direction), legs 1 and 3 seaward (Fig. 5.5)
    c, s = math.cos(math.radians(ang)), math.sin(math.radians(ang))
    fair_body = [r_fair * c, r_fair * s, round(z_fair_body, 5)]
    anchor = [cg_ref[0] + (r_fair + horiz) * c, cg_ref[1] + (r_fair + horiz) * s, 0.02]    # anchor plate 0.02 m high (Fig. 5.4)
    return fair_body, [round(v, 5) for v in anchor]

def config(k):
    cfg = {"_comment": f"Paredes (2016) moored-buoy benchmark, configuration {k.upper()}. Sources: Tables 5.3-5.14, Figs 5.4-5.5; ASSUMPTIONs listed in docs/paredes_validation.md. SI units.",
           "body": body(k),
           "environment": {"hydro": True, "seabed": k == "cat", "rho_w_kg_m3": rho, "water_surface_z_m": depth, "seabed_z_m": 0.02 if k == "cat" else 0.0},
           "soil": {"stiffness_Pa_per_m": 3.0e9, "damping_factor": 1.0, "friction": 0.3, "v_lim_m_s": 0.01},   # CAT bed: ASSUMPTION (values of the Chalmers chain)
           "lines": [], "initial": {"equilibrium": True},
           "numerics": {"dt_body_s": 0.002, "t_end_s": 40.0, "scheme": "rk4", "cfl": 0.5},
           "output": {"directory": f"out/paredes_{k}", "tag": k, "dt_out_s": 0.02, "write": False}}
    if k in ("con1", "con2"):
        cfg["line_defaults"] = {"length_m": 2.285, "segments": 23, "EA_N": 1.6e5, "mass_per_length_kg_m": 3.2e-3,        # Table 5.9
                                "weight_per_length_N_m": 8.0e-3, "hydro_diameter_m": 1.742e-3, "soil_diameter_m": 1.742e-3,   # D: equivalent-mass diameter, ASSUMPTION
                                "Cm": 3.8, "Cdt": 0.5, "Cdn": 2.5}                                                         # coefficients of Table 3.3 (chain): ASSUMPTION for the cable
        for i in range(3):
            fb, an = leg_geometry(i, 1.800)
            pts = [{"type": "floater", "arclength_m": 0.700, "mass_kg": floaters[i][0], "buoyancy_N": floaters[i][1],
                    "diameter_m": D_floater, "Cd": 0.5}]                                                                   # sphere Cd: ASSUMPTION
            if k == "con2":
                pts.append({"type": "clump", "arclength_m": 1.600, "mass_kg": clumps[i][0], "submerged_weight_N": clumps[i][1],
                            "diameter_m": 0.04, "Cd": 1.0})                                                                # clump size and Cd: ASSUMPTION
            cfg["lines"].append({"name": f"leg{i+1}", "anchor_m": an, "fairlead_body_m": fb, "point_elements": pts})
    elif k == "cat":
        cfg["line_defaults"] = {"length_m": 6.95, "segments": 40, "EA_N": 1.6e6, "mass_per_length_kg_m": 0.1447,           # Table 5.10; length: 5.0 m lying + 1.95 m suspended (Fig. 5.4)
                                "weight_per_length_N_m": 1.243, "hydro_diameter_m": 4.86e-3, "soil_diameter_m": 4.86e-3,    # equivalent-mass steel diameter, ASSUMPTION
                                "Cm": 3.8, "Cdt": 0.5, "Cdn": 2.5}
        for i in range(3):
            fb, an = leg_geometry(i, 6.660)
            cfg["lines"].append({"name": f"leg{i+1}", "anchor_m": an, "fairlead_body_m": fb, "initial": {"shape": "touchdown"}})
    else:                                                                                                                 # free buoy
        cfg["initial"] = {"equilibrium": False}
        cfg["numerics"] = {"dt_body_s": 0.001, "t_end_s": 20.0}
    return cfg

os.makedirs(out, exist_ok=True)
for k in ("free", "con1", "con2", "cat"):
    name = "free_buoy" if k == "free" else k
    json.dump(config(k), open(os.path.join(out, name + ".json"), "w"), indent=1)
    print("wrote", name, "C33 = %.2f N/m" % C33)
