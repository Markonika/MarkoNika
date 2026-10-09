"""Controller interface and the first two controllers.

A controller is sampled at its own period and holds its output between samples (zero-order
hold). It receives the measurement vector y = [phi, phi_dot, alpha, alpha_dot] and returns
the precession torque u [N m]. MPC, sliding mode and learning controllers (Layer 3) plug in
through the same `step` method.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np


class Controller(Protocol):
    def reset(self) -> None: ...
    def step(self, t: float, y: np.ndarray) -> float: ...


class Passive:
    """No actuation: the gyro precesses freely against its damper."""

    def reset(self) -> None:
        pass

    def step(self, t: float, y: np.ndarray) -> float:
        return 0.0


class SaturatedRateFeedback:
    """u = sat(k * phi_dot): torque on the precession axis proportional to platform pitch rate.

    Placeholder baseline, NOT a reproduction of any published controller. The sign and gain
    are to be tuned against the passive case; the saturation is the actuator torque limit.
    """

    def __init__(self, k: float, u_max: float) -> None:
        self.k = k
        self.u_max = u_max

    def reset(self) -> None:
        pass

    def step(self, t: float, y: np.ndarray) -> float:
        return float(np.clip(self.k * y[1], -self.u_max, self.u_max))
