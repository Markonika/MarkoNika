#!/usr/bin/env python3
"""Combined sensitivity study of the Paredes CON2 wave response to the uncertain inputs (docs/bem_paredes.md): pitch inertia I55, clump and floater node (position
along the rope), floater drag scale and line drag scale. Writes variant configuration directories and a job list; a SENSITIVITY STUDY, not a calibration (nothing is fitted).
usage: paredes_con2_sensitivity.py OUT_DIR [n_lhs=16] [seed=1]
Ranges: I55 in [0.87, 1.00] kg m^2 (thesis value without fixtures .. value implied by the measured pitch period with the BEM A55); clump node 15..17 (as run: 16);
floater node 6..8 (as run: 7); floater Cd scale 0..2; line Cdt/Cdn scale 0.5..2 (as run: 1)."""
import json, os, shutil, sys
import numpy as np

out = sys.argv[1]
nlhs = int(sys.argv[2]) if len(sys.argv) > 2 else 16
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "examples", "paredes")
base = json.load(open(os.path.join(src, "con2.json")))
os.makedirs(out, exist_ok=True)

def variant(name, I55, clump, floater, fl_drag, line_drag):
    d = os.path.join(out, "v_" + name)
    os.makedirs(d, exist_ok=True)
    for f in os.listdir(src): shutil.copy(os.path.join(src, f), d)
    c = json.loads(json.dumps(base))
    c["body"]["inertia_diag_kg_m2"] = [I55, I55, c["body"]["inertia_diag_kg_m2"][2]]
    c["line_defaults"]["Cdt"] *= line_drag; c["line_defaults"]["Cdn"] *= line_drag
    for l in c["lines"]:
        for p in l["point_elements"]:
            p.pop("arclength_m", None)
            if p["type"] == "clump": p["node"] = int(clump)
            else: p["node"] = int(floater); p["Cd"] *= fl_drag
    json.dump(c, open(os.path.join(d, "con2.json"), "w"), indent=1)

jobs = []   # (name, period)
def add(name, I55, clump, floater, fl_drag, line_drag, periods):
    variant(name, I55, clump, floater, fl_drag, line_drag)
    for T in periods: jobs.append((name, T))

# stage 1: corners of the box (lowest-pitch corner, its as-configured-drag neighbour, highest-pitch corner, as run)
add("asrun", 0.87, 16, 7, 1.0, 1.0, (1.30, 1.40))
add("cornerLow", 0.87, 15, 8, 0.0, 0.5, (1.30, 1.40))
add("cornerLowDrag1", 0.87, 15, 8, 1.0, 1.0, (1.30, 1.40))
add("cornerHigh", 1.00, 17, 6, 2.0, 2.0, (1.30, 1.40))
# stage 2: Latin hypercube over the box
rng = np.random.default_rng(seed)
def strat(n): return (rng.permutation(n) + rng.random(n)) / n
u = [strat(nlhs) for _ in range(5)]
for k in range(nlhs):
    I55 = 0.87 + 0.13 * u[0][k]
    clump = 15 + min(2, int(u[1][k] * 3)); floater = 6 + min(2, int(u[2][k] * 3))
    fl = 2.0 * u[3][k]; ln = 0.5 + 1.5 * u[4][k]
    add("lhs%02d" % k, round(I55, 4), clump, floater, round(fl, 3), round(ln, 3), (1.30,))
    print("lhs%02d  I55 %.3f clump node %d floater node %d floater drag x%.2f line drag x%.2f" % (k, I55, clump, floater, fl, ln))
with open(os.path.join(out, "jobs.txt"), "w") as f:
    for n, T in jobs: f.write("%s %.2f\n" % (n, T))
print(len(jobs), "jobs written to", os.path.join(out, "jobs.txt"))
