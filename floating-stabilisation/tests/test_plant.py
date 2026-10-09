import numpy as np

from fowt_stab import Passive, PitchGyroParams, energy, simulate, wave_moment_series
from fowt_stab.waves import jonswap


def zero(t):
    return 0.0


def test_energy_conserved_without_dissipation():
    # n=1, no damping, no torque: the gyroscopic coupling is skew-symmetric.
    p = PitchGyroParams(n=1, Cp=0.0, zeta_F=0.0)
    x0 = np.array([np.deg2rad(2.0), 0.0, 0.0, 0.0])
    res = simulate(p, Passive(), zero, t_end=200.0, dt=0.005, Ts=0.05, x0=x0, limit=False)
    E = np.array([energy(x, p) for x in res.x])
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-6


def test_free_decay_period_with_gyro_off():
    p = PitchGyroParams(rpm=0.0)                      # no spin: decoupled platform
    x0 = np.array([np.deg2rad(3.0), 0.0, 0.0, 0.0])
    res = simulate(p, Passive(), zero, t_end=200.0, dt=0.01, x0=x0)
    phi = res.x[:, 0]
    up = np.where((phi[:-1] < 0) & (phi[1:] >= 0))[0]            # upward zero crossings
    t_up = res.t[up] - phi[up] * (res.t[up + 1] - res.t[up]) / (phi[up + 1] - phi[up])
    period = np.mean(np.diff(t_up))
    expected = p.T_n / np.sqrt(1 - p.zeta_F**2)
    assert abs(period - expected) / expected < 1e-3


def test_rk4_fourth_order_convergence():
    p = PitchGyroParams(n=1, Cp=0.0, zeta_F=0.0)
    x0 = np.array([np.deg2rad(2.0), 0.0, 0.0, 0.0])

    def final(dt):
        return simulate(p, Passive(), zero, t_end=20.0, dt=dt, Ts=0.1, x0=x0, limit=False).x[-1]

    ref = final(0.0025)
    e1 = np.linalg.norm(final(0.02) - ref)
    e2 = np.linalg.norm(final(0.01) - ref)
    assert 12.0 < e1 / e2 < 20.0                       # ~16 for a 4th-order scheme


def test_precession_limit_respected():
    p = PitchGyroParams(Cp=0.0, alpha_max=np.deg2rad(20.0))
    t = np.arange(0, 600.0, 0.01)
    M = wave_moment_series(t, Hs=4.0, Tp=12.0, seed=1)
    res = simulate(p, Passive(), lambda s: np.interp(s, t, M), t_end=500.0)
    assert np.max(np.abs(res.x[:, 2])) <= p.alpha_max + 1e-12


def test_jonswap_hits_target_hs():
    w = np.linspace(0.2, 2.5, 4000)
    S = jonswap(w, Hs=4.0, Tp=12.0)
    assert abs(4.0 * np.sqrt(np.trapezoid(S, w)) - 4.0) < 1e-6


def test_passive_gyro_reduces_pitch_in_irregular_waves():
    t = np.arange(0, 1800.0, 0.01)
    M = wave_moment_series(t, Hs=4.0, Tp=12.0, seed=3)
    moment = lambda s: np.interp(s, t, M)
    on = simulate(PitchGyroParams(), Passive(), moment, t_end=1700.0)
    off = simulate(PitchGyroParams(rpm=0.0), Passive(), moment, t_end=1700.0)
    keep = on.t > 200.0                                # drop the transient
    assert on.x[keep, 0].std() < off.x[keep, 0].std()
