#!/usr/bin/env python3
"""Combined sensitivity study of the Paredes moored-buoy wave response (docs/bem_paredes.md) to the uncertain inputs. Writes variant configuration directories and a job
list (lines "<variant> <period>"); a SENSITIVITY STUDY, not a calibration (nothing is fitted).
usage: paredes_sensitivity.py OUT_DIR [config=con2|con1|cat] [n_lhs=16] [seed=1]
Common: pitch inertia I55 in [0.87, 1.00] kg m^2 (thesis value without fixtures .. value implied by the measured pitch period with the BEM A55); line Cdt/Cdn scale 0.5..2.
con2: clump node 15..17 (as run 16), floater node 6..8 (as run 7), floater Cd scale 0..2 (+ four corner cases, periods 1.3 and 1.4 s for the corners, 1.3 s for the sample).
con1: floater node 6..8, floater Cd scale 0..2 (both periods for the sample).
cat : total chain length 6.95..7.00 m (Fig. 5.4 labels give 6.95 m, Table 5.6 implies 7.0 m); no point elements (both periods for the sample)."""
import json, os, shutil, sys
import numpy as np

out = sys.argv[1]
name = sys.argv[2] if len(sys.argv) > 2 else "con2"
nlhs = int(sys.argv[3]) if len(sys.argv) > 3 else 16
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "examples", "paredes")
base = json.load(open(os.path.join(src, name + ".json")))
os.makedirs(out, exist_ok=True)

def variant(vn, I55, clump=None, floater=None, fl_drag=1.0, line_drag=1.0, length=None):
    d = os.path.join(out, "v_" + vn)
    os.makedirs(d, exist_ok=True)
    for f in os.listdir(src): shutil.copy(os.path.join(src, f), d)
    c = json.loads(json.dumps(base))
    c["body"]["inertia_diag_kg_m2"] = [I55, I55, c["body"]["inertia_diag_kg_m2"][2]]
    c["line_defaults"]["Cdt"] *= line_drag; c["line_defaults"]["Cdn"] *= line_drag
    if length is not None: c["line_defaults"]["length_m"] = length
    for l in c["lines"]:
        for p in l.get("point_elements", []):
            p.pop("arclength_m", None)
            if p["type"] == "clump": p["node"] = int(clump)
            else: p["node"] = int(floater); p["Cd"] *= fl_drag
    json.dump(c, open(os.path.join(d, name + ".json"), "w"), indent=1)

jobs = []
def add(vn, periods, **kw):
    variant(vn, **kw)
    for T in periods: jobs.append((vn, T))

rng = np.random.default_rng(seed)
def strat(n): return (rng.permutation(n) + rng.random(n)) / n
nu = {"con2": 5, "con1": 4, "cat": 3}[name]
u = [strat(nlhs) for _ in range(nu)]
if name == "con2":
    add("asrun", (1.30, 1.40), I55=0.87, clump=16, floater=7)
    add("cornerLow", (1.30, 1.40), I55=0.87, clump=15, floater=8, fl_drag=0.0, line_drag=0.5)
    add("cornerLowDrag1", (1.30, 1.40), I55=0.87, clump=15, floater=8)
    add("cornerHigh", (1.30, 1.40), I55=1.00, clump=17, floater=6, fl_drag=2.0, line_drag=2.0)
    per = (1.30,)
elif name == "con1":
    add("asrun", (1.30, 1.40), I55=0.87, floater=7)
    per = (1.30, 1.40)
else:
    add("asrun", (1.30, 1.40), I55=0.87)
    per = (1.30, 1.40)
for k in range(nlhs):
    I55 = round(0.87 + 0.13 * u[0][k], 4)
    if name == "con2":
        clump = 15 + min(2, int(u[1][k] * 3)); floater = 6 + min(2, int(u[2][k] * 3)); fl = round(2.0 * u[3][k], 3); ln = round(0.5 + 1.5 * u[4][k], 3)
        add("lhs%02d" % k, per, I55=I55, clump=clump, floater=floater, fl_drag=fl, line_drag=ln)
        print("lhs%02d  I55 %.3f clump node %d floater node %d floater drag x%.2f line drag x%.2f" % (k, I55, clump, floater, fl, ln))
    elif name == "con1":
        floater = 6 + min(2, int(u[1][k] * 3)); fl = round(2.0 * u[2][k], 3); ln = round(0.5 + 1.5 * u[3][k], 3)
        add("lhs%02d" % k, per, I55=I55, floater=floater, fl_drag=fl, line_drag=ln)
        print("lhs%02d  I55 %.3f floater node %d floater drag x%.2f line drag x%.2f" % (k, I55, floater, fl, ln))
    else:
        ln = round(0.5 + 1.5 * u[1][k], 3); L = round(6.95 + 0.05 * u[2][k], 4)
        add("lhs%02d" % k, per, I55=I55, line_drag=ln, length=L)
        print("lhs%02d  I55 %.3f line drag x%.2f length %.4f m" % (k, I55, ln, L))
with open(os.path.join(out, "jobs.txt"), "w") as f:
    for n, T in jobs: f.write("%s %.2f\n" % (n, T))
print(len(jobs), "jobs written to", os.path.join(out, "jobs.txt"))
