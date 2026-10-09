#!/usr/bin/env python3
"""Lopez-Olocco et al. (2022) comparison: max/min fairlead tension for amplitudes x periods of one configuration.
usage: lopez_compare.py OUTDIR CONFIG_JSON CONFIGURATION(WO_CW|CW1|CW2) [key=value ...]
Reads examples/lopezolocco/measured.csv (Tables 8, 9 of the paper); writes OUTDIR/results.json and prints the comparison."""
import csv, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
exe = os.path.join(root, "build", "mooring_run")
out, cfgfile, conf = sys.argv[1:4]; extra = sys.argv[4:]
os.makedirs(out, exist_ok=True)
all_periods = [2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]; all_amps = [0.125, 0.15, 0.175, 0.2, 0.225]
periods = [float(x) for x in os.environ["LO_PERIODS"].split(",")] if "LO_PERIODS" in os.environ else all_periods   # optional subset
amps = [float(x) for x in os.environ["LO_AMPS"].split(",")] if "LO_AMPS" in os.environ else all_amps
meas = {}
for r in csv.reader(l for l in open(os.path.join(root, "examples", "lopezolocco", "measured.csv")) if not l.startswith("#")):
    if r[0] in ("max", "min") and r[1] == conf:
        for p, v in zip(all_periods, r[3:]): meas[(r[0], float(r[2]), p)] = float(v)

def run(key):
    A, T = key; tag = f"{conf}_A{A}_T{T}"
    subprocess.run([exe, cfgfile, f"motion.radius_m={A}", f"motion.period_s={T}", f"output.directory={out}/runs", f"output.tag={tag}"] + extra, capture_output=True, text=True)
    rows = list(csv.DictReader(open(f"{out}/runs/{tag}_timeseries.csv")))
    cyc = [float(l.split(",")[1]) for l in open(f"{out}/runs/{tag}_cycles.csv").read().split("\n")[1:] if l]
    n = len(cyc) - 1
    seg = [float(r["top_tension_N"]) for r in rows if (n - 4) * T <= float(r["time_s"]) <= n * T]
    return key, {"max": sum(cyc[n - 4:n]) / 4.0, "min": min(seg)}

with ThreadPoolExecutor(max_workers=4) as ex: res = dict(ex.map(run, [(a, t) for a in amps for t in periods]))
json.dump({f"{a}_{t}": v for (a, t), v in res.items()}, open(f"{out}/results.json", "w"), indent=1)
dm = []; dn = []
print(f"{conf}: model / measured, differences in % (max) and N (min)")
print("   A \\ T " + "".join(f"{t:>14}" for t in periods))
for a in amps:
    print(f"{a:7.3f} " + "".join(f"  {res[(a,t)]['max']:5.2f}/{meas[('max',a,t)]:5.2f}" for t in periods) + "   max")
    print(f"{'':7} " + "".join(f"  {res[(a,t)]['min']:5.2f}/{meas[('min',a,t)]:5.2f}" for t in periods) + "   min")
    for t in periods:
        dm.append(100 * (res[(a, t)]["max"] - meas[("max", a, t)]) / meas[("max", a, t)]); dn.append(res[(a, t)]["min"] - meas[("min", a, t)])
print(f"max tension: mean diff {sum(dm)/len(dm):+.2f} %, mean |diff| {sum(abs(x) for x in dm)/len(dm):.2f} %, worst {max(dm, key=abs):+.2f} %")
print(f"min tension: mean diff {sum(dn)/len(dn):+.3f} N, mean |diff| {sum(abs(x) for x in dn)/len(dn):.3f} N, worst {max(dn, key=abs):+.3f} N")
