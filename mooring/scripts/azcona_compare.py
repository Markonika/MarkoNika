#!/usr/bin/env python3
"""Run the Azcona et al. (2017) cases (config 1/2 x T = 1.58, 3.16, 4.74 s) and the quasi-static curves with mooring_run.
usage: azcona_compare.py OUTDIR [key=value ...]  (overrides go to every run). Writes OUTDIR/azcona_results.json and prints a table."""
import csv, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
exe = os.path.join(root, "build", "mooring_run")
out = sys.argv[1]; extra = sys.argv[2:]
os.makedirs(out, exist_ok=True)
# measured, read by eye from Figs 5/6 of the paper (mean markers; +-0.5 N); static from Table 5
meas = {("conf1", 1.58): (14.5, 3.0), ("conf1", 3.16): (9.8, 6.1), ("conf1", 4.74): (9.2, 7.3),
        ("conf2", 1.58): (44.5, 1.5), ("conf2", 3.16): (23.5, 7.5), ("conf2", 4.74): (19.5, 11.0)}
static_meas = {"conf1": 8.13, "conf2": 14.48}
d = {"conf1": 19.364, "conf2": 19.872}

def run(key):
    conf, T = key; tag = f"{conf}_T{T}"
    cfg = os.path.join(root, "examples", "azcona", conf + ".json")
    subprocess.run([exe, cfg, f"motion.period_s={T}", f"output.directory={out}/runs", f"output.tag={tag}"] + extra, capture_output=True, text=True)
    rows = list(csv.DictReader(open(f"{out}/runs/{tag}_timeseries.csv")))
    t = [float(r["time_s"]) for r in rows]; ten = [float(r["top_tension_N"]) for r in rows]
    cyc = [float(l.split(",")[1]) for l in open(f"{out}/runs/{tag}_cycles.csv").read().split("\n")[1:] if l]
    ncyc = len(cyc) - 1
    t0 = (ncyc - 4) * T                                   # last 4 complete cycles
    seg = [x for tt, x in zip(t, ten) if tt >= t0 and tt <= ncyc * T]
    mx = sum(cyc[ncyc - 4:ncyc]) / 4.0                     # per-step maxima of the last 4 cycles (the CSV is decimated)
    der = json.load(open(f"{out}/runs/{tag}_params.json"))["derived"]
    return key, {"max": mx, "min": min(seg), "mean": sum(seg) / len(seg), "dt": der["dt_used_s"], "slack": der["slack_segment_evals"]}

def static(conf, x):
    cfg = os.path.join(root, "examples", "azcona", conf + ".json")
    p = subprocess.run([exe, cfg, "motion.type=none", "output.write=false", f"fairlead_rest_m=[{x},0,5.0]"] + [e for e in extra if not e.startswith("motion")], capture_output=True, text=True)
    return float(p.stdout.split()[3])

with ThreadPoolExecutor(max_workers=4) as ex: res = dict(ex.map(run, sorted(meas)))
amp = 0.125
for e in extra:
    if e.startswith("motion.radius_m="): amp = float(e.split("=")[1])
res_s = {c: {"static": static(c, d[c]), "qs_lo": static(c, d[c] - amp), "qs_hi": static(c, d[c] + amp)} for c in d}
json.dump({"dynamic": {f"{k[0]}_{k[1]}": v for k, v in res.items()}, "static": res_s, "overrides": extra}, open(f"{out}/azcona_results.json", "w"), indent=1)
print("static top tension: " + ", ".join(f"{c}: {res_s[c]['static']:.2f} N (measured {static_meas[c]})" for c in d))
print("quasi-static ends: " + ", ".join(f"{c}: {res_s[c]['qs_lo']:.2f} .. {res_s[c]['qs_hi']:.2f} N" for c in d))
print(f"{'case':>12} {'T':>5} | {'max model':>9} {'max meas':>8} {'diff':>7} | {'min model':>9} {'min meas':>8}")
for (c, T), (mm, mn) in sorted(meas.items()):
    r = res[(c, T)]
    print(f"{c:>12} {T:5.2f} | {r['max']:9.2f} {mm:8.1f} {100*(r['max']-mm)/mm:+6.1f}% | {r['min']:9.2f} {mn:8.1f}")
