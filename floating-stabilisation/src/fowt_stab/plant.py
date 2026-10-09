"""Layer-1 reduced plant: platform pitch coupled to a counter-rotating gyro array.

State x = [phi, phi_dot, alpha, alpha_dot]
    phi    platform pitch angle [rad]
    alpha  gyro precession angle [rad]

Equations (n identical gyros lumped into one equivalent rotor; the counter-rotating
pair cancels the yaw moment, so only the pitch/precession pair is kept):

    (I_F + n*I_g) phi'' + C_F phi' + K_F phi = M_wave - n*Izz*w_g*alpha'*cos(alpha)
    n*Ixx alpha'' + n*Cp alpha'              = n*Izz*w_g*phi'*cos(alpha) + u

This is the form used by Wang C. et al. (Ocean Eng. 329, 121147, Eq. 5), with the
small-angle-in-phi assumption but the full cos(alpha) kept. For n = 1, Cp = 0, C_F = 0,
u = 0 the gyroscopic coupling is skew-symmetric, so total energy is conserved. That is
the first verification target (tests/test_plant.py).

Parameters below are PLACEHOLDERS shaped like the OC4-DeepCwind case in the papers
(pitch inertia 8.35e9 kg m^2, ~25 s free-decay period, gyro Izz 1.5e5 kg m^2). They are
not a calibrated model; the platform stiffness is set only to give the 25 s period.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np


@dataclass(frozen=True)
class PitchGyroParams:
    # platform (effective pitch inertia incl. added mass, placeholder)
    I_F: float = 8.35e9          # kg m^2
    T_n: float = 25.0            # s, target pitch natural period (gyro off)
    zeta_F: float = 0.02         # -, hydrodynamic damping ratio, placeholder
    # gyro (per unit)
    Izz: float = 1.5e5           # kg m^2, spin inertia
    Ixx: float = 7.5e4           # kg m^2, precession-axis inertia (placeholder ~ Izz/2)
    rpm: float = 6000.0          # spin speed
    Cp: float = 8.0e5            # N m s/rad, precession damping (Wang C. default)
    n: int = 2                   # number of units (counter-rotating pair)
    alpha_max: float = np.deg2rad(70.0)   # precession angle limit

    @property
    def K_F(self) -> float:
        return self.I_F * (2.0 * np.pi / self.T_n) ** 2

    @property
    def C_F(self) -> float:
        return 2.0 * self.zeta_F * np.sqrt(self.K_F * self.I_F)

    @property
    def w_g(self) -> float:
        return self.rpm * 2.0 * np.pi / 60.0

    def with_(self, **kw) -> "PitchGyroParams":
        return replace(self, **kw)


def rhs(x: np.ndarray, M_wave: float, u: float, p: PitchGyroParams) -> np.ndarray:
    """Time derivative of the state for given wave moment [N m] and precession torque u [N m]."""
    phi, phid, alpha, alphad = x
    h = p.n * p.Izz * p.w_g                       # total spin angular momentum
    ca = np.cos(alpha)
    phidd = (M_wave - p.C_F * phid - p.K_F * phi - h * alphad * ca) / (p.I_F + 0.0)
    alphadd = (h * phid * ca + u - p.n * p.Cp * alphad) / (p.n * p.Ixx)
    return np.array([phid, phidd, alphad, alphadd])


def energy(x: np.ndarray, p: PitchGyroParams) -> float:
    """Mechanical energy of the conservative part (used for verification)."""
    phi, phid, _, alphad = x
    return 0.5 * p.I_F * phid**2 + 0.5 * p.K_F * phi**2 + 0.5 * p.n * p.Ixx * alphad**2


def apply_precession_limit(x: np.ndarray, p: PitchGyroParams) -> np.ndarray:
    """Hard stop on the precession angle: clamp and remove outward velocity.

    Applied after each integrator step (non-smooth, so it is kept outside the RK stages).
    """
    x = x.copy()
    if x[2] > p.alpha_max:
        x[2] = p.alpha_max
        x[3] = min(x[3], 0.0)
    elif x[2] < -p.alpha_max:
        x[2] = -p.alpha_max
        x[3] = max(x[3], 0.0)
    return x
