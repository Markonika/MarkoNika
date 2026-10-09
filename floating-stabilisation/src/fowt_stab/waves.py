"""Irregular-wave excitation for the Layer-1 plant.

JONSWAP spectrum, random-phase superposition. Frequencies are spaced evenly and then
jittered inside each bin so the record does not repeat after 2*pi/d_omega (a known flaw
of equally spaced components).

The pitch-moment transfer function here is a PLACEHOLDER (flat times a roll-off), to be
replaced by the BEM-derived excitation RAO (Capytaine/NEMOH) in Layer 2.
"""

from __future__ import annotations

import numpy as np


def jonswap(omega: np.ndarray, Hs: float, Tp: float, gamma: float = 3.3) -> np.ndarray:
    """JONSWAP spectral density S(omega) [m^2 s/rad], normalised to the target Hs."""
    wp = 2.0 * np.pi / Tp
    sigma = np.where(omega <= wp, 0.07, 0.09)
    r = np.exp(-((omega - wp) ** 2) / (2.0 * sigma**2 * wp**2))
    S = omega**-5 * np.exp(-1.25 * (wp / omega) ** 4) * gamma**r
    m0 = np.trapezoid(S, omega)
    return S * (Hs**2 / 16.0) / m0


def wave_moment_series(
    t: np.ndarray,
    Hs: float,
    Tp: float,
    seed: int,
    moment_rao: float = 6.0e8,
    w_roll: float = 1.2,
    n_comp: int = 200,
    w_min: float = 0.2,
    w_max: float = 2.5,
) -> np.ndarray:
    """Wave pitch moment M(t) [N m] for the time vector t [s]."""
    rng = np.random.default_rng(seed)
    edges = np.linspace(w_min, w_max, n_comp + 1)
    w = edges[:-1] + rng.uniform(0.0, 1.0, n_comp) * np.diff(edges)   # jittered
    dw = np.diff(edges)
    S = jonswap(w, Hs, Tp)
    amp = np.sqrt(2.0 * S * dw)                                        # elevation amplitude
    eps = rng.uniform(0.0, 2.0 * np.pi, n_comp)
    H = moment_rao * np.exp(-((w / w_roll) ** 2))                     # N m per m of wave
    return (amp * H) @ np.cos(np.outer(w, t) + eps[:, None])
