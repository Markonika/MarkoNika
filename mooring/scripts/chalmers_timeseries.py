#!/usr/bin/env python3
"""Plot raw top/anchor tension for r=0.2 m, T=3.5 s and T=1.25 s (last 2 cycles, unfiltered).
usage: chalmers_timeseries.py OUTDIR"""
import csv, subprocess, sys, os
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
exe = os.path.join(root, "build", "mooring_run"); cfg = os.path.join(root, "examples", "chalmers", "chalmers_config.json")
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
for a, T in zip(ax, (3.5, 1.25)):
    tag = f"ts_r0.2_T{T}"
    subprocess.run([exe, cfg, "motion.radius_m=0.2", f"motion.period_s={T}", f"output.directory={out}", f"output.tag={tag}", "output.dt_out_s=0.001", "motion.cycles=10"], capture_output=True)
    rows = list(csv.DictReader(open(f"{out}/{tag}_timeseries.csv")))
    t = [float(r["time_s"]) for r in rows]; top = [float(r["top_tension_N"]) for r in rows]; anc = [float(r["anchor_tension_N"]) for r in rows]
    i0 = next(i for i, x in enumerate(t) if x >= 8 * T)
    a.plot([(x - 8 * T) / T for x in t[i0:]], top[i0:], lw=0.8, label="top"); a.plot([(x - 8 * T) / T for x in t[i0:]], anc[i0:], lw=0.8, label="anchor")
    a.set_title(f"r = 0.2 m, T = {T} s (raw, unfiltered)"); a.set_xlabel("t / T (cycles 8-9)"); a.set_ylabel("tension [N]"); a.legend(fontsize=8)
    seg = top[i0:]; print(f"T={T}: min {min(seg):.2f} max {max(seg):.2f} N over cycles 8-9")
fig.tight_layout(); fig.savefig(f"{out}/timeseries.png", dpi=150)
