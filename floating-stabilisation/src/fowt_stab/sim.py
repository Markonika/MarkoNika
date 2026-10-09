"""Fixed-step RK4 simulation with a zero-order-hold controller.

The plant step dt is a divisor of the controller period Ts. The precession hard stop is applied
after each RK4 step (outside the stages).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .control import Controller
from .plant import PitchGyroParams, apply_precession_limit, rhs


@dataclass
class Result:
    t: np.ndarray
    x: np.ndarray        # (N, 4)
    u: np.ndarray        # (N,)


def rk4_step(x, M0, M_half, M1, u, dt, p):
    """One RK4 step with the wave moment sampled at t, t+dt/2, t+dt."""
    k1 = rhs(x, M0, u, p)
    k2 = rhs(x + 0.5 * dt * k1, M_half, u, p)
    k3 = rhs(x + 0.5 * dt * k2, M_half, u, p)
    k4 = rhs(x + dt * k3, M1, u, p)
    return x + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


def simulate(
    p: PitchGyroParams,
    controller: Controller,
    moment,                 # callable t -> N m
    t_end: float,
    dt: float = 0.01,
    Ts: float = 0.05,
    x0: np.ndarray | None = None,
    limit: bool = True,
) -> Result:
    n_sub = int(round(Ts / dt))
    if abs(n_sub * dt - Ts) > 1e-12:
        raise ValueError("Ts must be an integer multiple of dt")
    n_steps = int(round(t_end / dt))
    x = np.zeros(4) if x0 is None else np.asarray(x0, float).copy()
    t = np.arange(n_steps + 1) * dt
    X = np.empty((n_steps + 1, 4))
    U = np.empty(n_steps + 1)
    X[0] = x
    controller.reset()
    u = 0.0
    for i in range(n_steps):
        if i % n_sub == 0:
            u = controller.step(t[i], x)
        U[i] = u
        x = rk4_step(x, moment(t[i]), moment(t[i] + 0.5 * dt), moment(t[i] + dt), u, dt, p)
        if limit:
            x = apply_precession_limit(x, p)
        X[i + 1] = x
    U[-1] = u
    return Result(t, X, U)
