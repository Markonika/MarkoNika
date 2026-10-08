#!/usr/bin/env python3
"""Convergence of the Chalmers mean-maximum top tension in N (segments) and dt, plus RK4 vs Verlet.
usage: chalmers_convergence.py OUTDIR      -> OUTDIR/convergence.json, convergence.png"""
import json, subprocess, sys, os
from concurrent.futures import ThreadPoolExecutor
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
exe = os.path.join(root, "build", "mooring_run"); cfg = os.path.join(root, "examples", "chalmers", "chalmers_config.json")
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
cases = [(0.2, 1.25), (0.2, 3.5), (0.1, 2.0)]
Ns = [16, 33, 66, 132, 264]
dts = [1.1e-4, 5.5e-5, 2.75e-5, 1.4e-5]

def run(job):
    tag, rad, T, over = job
    a = [exe, cfg, f"motion.radius_m={rad}", f"motion.period_s={T}", f"output.directory={out}/runs", f"output.tag={tag}", "output.write=true", "output.dt_out_s=0.05"] + over
    subprocess.run(a, capture_output=True, text=True)
    d = json.load(open(f"{out}/runs/{tag}_params.json"))["derived"]
    return tag, {"mean_max_N": d["mean_cycle_max_N"], "dt_s": d["dt_used_s"], "wall_s": d["wall_s"], "finite": d["finite"], "static_top_N": d["static_top_tension_N"]}

jobs = []
for rad, T in cases:
    for N in Ns: jobs.append((f"N{N}_r{rad}_T{T}", rad, T, [f"line.segments={N}"]))
    for dt in dts: jobs.append((f"dt{dt:g}_r{rad}_T{T}", rad, T, [f"numerics.dt_s={dt}"]))
    jobs.append((f"verlet_r{rad}_T{T}", rad, T, ["numerics.scheme=verlet"]))
with ThreadPoolExecutor(max_workers=4) as ex: res = dict(ex.map(run, jobs))
json.dump(res, open(f"{out}/convergence.json", "w"), indent=1)
for rad, T in cases:
    print(f"\nr={rad} T={T}\n  N-study (dt = CFL/soil limit):")
    for N in Ns: r = res[f"N{N}_r{rad}_T{T}"]; print(f"   N={N:4d} l0={33/N:.3f} dt={r['dt_s']:.2e} max={r['mean_max_N']:.4f} static={r['static_top_N']:.4f} wall={r['wall_s']:.0f}s")
    print("  dt-study (N=66):")
    for dt in dts: r = res[f"dt{dt:g}_r{rad}_T{T}"]; print(f"   dt={dt:.2e} max={r['mean_max_N']:.4f}")
    v = res[f"verlet_r{rad}_T{T}"]; print(f"  Verlet (N=66, default dt): {v['mean_max_N']:.4f}")
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    for rad, T in cases:
        ax[0].plot([33/N for N in Ns], [res[f"N{N}_r{rad}_T{T}"]["mean_max_N"] for N in Ns], "o-", label=f"r={rad}, T={T}")
        ax[1].plot(dts, [res[f"dt{dt:g}_r{rad}_T{T}"]["mean_max_N"] for dt in dts], "o-", label=f"r={rad}, T={T}")
    ax[0].set_xlabel("segment length l0 [m]"); ax[0].set_ylabel("mean max top tension [N]"); ax[0].invert_xaxis(); ax[0].legend(fontsize=8)
    ax[1].set_xlabel("time step [s] (N=66)"); ax[1].set_xscale("log"); ax[1].invert_xaxis()
    fig.tight_layout(); fig.savefig(f"{out}/convergence.png", dpi=150)
except Exception as e: print("plot skipped:", e)
