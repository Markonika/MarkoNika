"""Passive-gyro sweep over precession damping and spin speed (Layer-1 plausibility check).

Compares only the TRENDS with Wang C. et al. (Ocean Eng. 329, 121147): reduction should grow
with spin speed with diminishing returns, and should have an interior optimum in Cp (poor at
very low and very high damping). Absolute numbers depend on the placeholder wave-moment model.

Run:  PYTHONPATH=src python3 examples/sweep_passive.py
"""

import numpy as np

from fowt_stab import Passive, PitchGyroParams, simulate, wave_moment_series

T_END, DT = 1800.0, 0.01
t = np.arange(0, T_END + 1, DT)
SEEDS = (1, 2, 3)
moments = [wave_moment_series(t, Hs=4.0, Tp=12.0, seed=s) for s in SEEDS]


def pitch_std(p, scale=1.0):
    out = []
    for M in moments:
        res = simulate(p, Passive(), lambda s, M=M: scale * np.interp(s, t, M), t_end=T_END, dt=DT)
        out.append(res.x[res.t > 300.0, 0].std())
    return np.mean(out)


# Calibrate the placeholder wave moment so the gyro-off response matches the ~0.5 deg pitch
# standard deviation seen in the papers. The gyro-off platform is linear, so one scale works.
TARGET = np.deg2rad(0.5)
SCALE = TARGET / pitch_std(PitchGyroParams(rpm=0.0))
pitch_std_ = pitch_std
pitch_std = lambda p: pitch_std_(p, SCALE)     # noqa: E731

base = pitch_std(PitchGyroParams(rpm=0.0))
print(f"wave-moment scale factor {SCALE:.4f} (calibration only)")
print(f"gyro off: pitch std = {np.rad2deg(base):.3f} deg (mean of {len(SEEDS)} seeds)")

print("\nspin sweep (Cp = 8e5 N m s/rad)")
for rpm in (2000, 4000, 6000, 8000, 12000):
    s = pitch_std(PitchGyroParams(rpm=rpm))
    print(f"  {rpm:6d} rpm: reduction {100 * (1 - s / base):5.1f} %")

print("\nCp sweep (6000 rpm)")
for Cp in (4e3, 4e4, 2e5, 8e5, 1.4e6):
    s = pitch_std(PitchGyroParams(Cp=Cp))
    print(f"  Cp={Cp:8.1e}: reduction {100 * (1 - s / base):5.1f} %")
