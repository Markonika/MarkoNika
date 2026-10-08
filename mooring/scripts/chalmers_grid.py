#!/usr/bin/env python3
"""Run the 30-case Chalmers grid with mooring_run and compare with Table 7.

usage: chalmers_grid.py OUTDIR [key=value ...]     (overrides are passed to every case)
Writes OUTDIR/grid_results.csv and OUTDIR/grid_regression.json (+ PNG if matplotlib is present).
Analysis script only; the solver does not depend on Python.
"""
import csv, json, subprocess, sys, os
from concurrent.futures import ThreadPoolExecutor

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
exe = os.path.join(root, "build", "mooring_run")
cfg = os.path.join(root, "examples", "chalmers", "chalmers_config.json")
table = os.path.join(root, "examples", "chalmers", "chalmers_table7_max_tension.csv")
out = sys.argv[1]; extra = sys.argv[2:]
os.makedirs(out, exist_ok=True)

rows = [l for l in open(table) if l.strip() and not l.startswith("#")]
rd = list(csv.reader(rows)); hdr = rd[0]
radii = [float(h.split("_")[1]) for h in hdr[1:]]
meas = {}
for r in rd[1:]:
    for rad, v in zip(radii, r[1:]):
        meas[(float(r[0]), rad)] = float(v)

def run(key):
    T, rad = key
    tag = f"r{rad:.3f}_T{T:.2f}"
    args = [exe, cfg, f"motion.radius_m={rad}", f"motion.period_s={T}", f"output.directory={out}/runs",
            f"output.tag={tag}", "output.dt_out_s=0.02"] + extra
    p = subprocess.run(args, capture_output=True, text=True)
    pj = json.load(open(f"{out}/runs/{tag}_params.json"))["derived"]
    return key, pj, p.returncode

with ThreadPoolExecutor(max_workers=4) as ex:
    res = list(ex.map(run, sorted(meas)))

n = len(res); xs = []; ys = []
with open(f"{out}/grid_results.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["period_s", "radius_m", "measured_N", "simulated_N", "rel_diff", "cycles_averaged", "slack_evals", "dt_s", "finite"])
    for (T, rad), pj, rc in res:
        m = meas[(T, rad)]; s = pj["mean_cycle_max_N"]
        xs.append(m); ys.append(s)
        w.writerow([T, rad, m, f"{s:.4f}", f"{(s-m)/m:.4f}", pj["cycles_averaged"], pj["slack_segment_evals"], pj["dt_used_s"], pj["finite"]])

mx = sum(xs)/n; my = sum(ys)/n
sxx = sum((x-mx)**2 for x in xs); syy = sum((y-my)**2 for y in ys); sxy = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
slope = sxy/sxx; icpt = my - slope*mx; r2 = sxy**2/(sxx*syy)
# r^2 of the 1:1 line (coefficient of determination of simulated vs measured, no fit)
ss_res = sum((y-x)**2 for x, y in zip(xs, ys)); r2_11 = 1 - ss_res/sxx
rms = (ss_res/n)**0.5; bias = sum(y-x for x, y in zip(xs, ys))/n
summ = {"n": n, "r2_regression": r2, "slope": slope, "intercept": icpt, "r2_vs_1to1": r2_11,
        "rmse_N": rms, "mean_bias_N": bias, "mean_rel_diff": sum((y-x)/x for x, y in zip(xs, ys))/n,
        "overrides": extra}
json.dump(summ, open(f"{out}/grid_regression.json", "w"), indent=2)
print(json.dumps(summ, indent=2))
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.2, 5))
    ax.errorbar(xs, ys, xerr=[0.05*x for x in xs], fmt="o", ms=4, capsize=2, label="cases (±5 % reading error)")
    lim = [min(xs+ys)*0.9, max(xs+ys)*1.05]; ax.plot(lim, lim, "k--", lw=1, label="1:1")
    ax.plot(lim, [icpt+slope*v for v in lim], "r-", lw=1, label=f"fit: y={slope:.2f}x{icpt:+.1f}, r²={r2:.3f}")
    ax.set_xlabel("measured mean max tension [N]"); ax.set_ylabel("simulated [N]"); ax.legend(fontsize=8); ax.set_aspect("equal")
    fig.tight_layout(); fig.savefig(f"{out}/grid_scatter.png", dpi=150)
except Exception as e:
    print("plot skipped:", e)
