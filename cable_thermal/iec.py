"""Standard methods used as independent cross-checks.

IEC 60287-1-1 / 60287-2-1 (steady-state rating of a single cable in free air)
--------------------------------------------------------------------------------
For a single-core LV cable without metallic sheath/armour (lambda1 = lambda2 = 0,
no dielectric losses, no T2):

    I = sqrt( dtheta / ( R (T1 + T3 + T4) ) )

    T1 = rho_T,ins / (2 pi) ln(1 + 2 t1 / dc)              insulation
    T3 = rho_T,jac / (2 pi) ln(1 + 2 t3 / Da)              jacket (Da = diameter under jacket)
    T4 = 1 / (pi De h (dtheta_s)^(1/4))                    free air, no sun
    h  = Z / De^g + E   (single cable: Z = 0.21, E = 3.94, g = 0.60)

    (dtheta_s)^(1/4)_{n+1} = [ dtheta / (1 + K_A (dtheta_s)^(1/4)_n) ]^(1/4),
    K_A = pi De h (T1 + T3),   initial (dtheta_s)^(1/4) = 2

derivation of the iteration: W_tot = pi De h dtheta_s^(5/4) (from T4) and
dtheta - dtheta_s = W_tot (T1 + T3)  =>  dtheta = dtheta_s (1 + K_A dtheta_s^(1/4)).

IEC 60949 (adiabatic short-circuit temperature rise)
----------------------------------------------------
    I^2 t = K^2 S^2 ln((theta_f + beta) / (theta_i + beta)),
    K = sqrt( sigma_c (beta + 20) / rho_20 )   [A s^0.5 / mm^2 with S in mm^2]

For copper: sigma_c = 3.45e6 J/(K m^3), beta = 234.5 K, rho_20 = 1.7241e-8 Ohm m -> K = 226.
The linear resistance law of the paper, r = r20 (1 + alpha (T - 20)), corresponds to
beta = 1/alpha - 20.
"""
from __future__ import annotations

import numpy as np

from .cables import CableSpec

IEC_FREE_AIR_SINGLE = dict(Z=0.21, E=3.94, g=0.60)


def iec60287_T1_T3(cable: CableSpec, k_ins: float, k_jac: float):
    T1 = 1.0 / k_ins / (2 * np.pi) * np.log(1 + 2 * cable.insulation_thickness / cable.conductor_diameter)
    T3 = 1.0 / k_jac / (2 * np.pi) * np.log(1 + 2 * cable.jacket_thickness / (2 * cable.R_B))
    return T1, T3


def iec60287_free_air(cable: CableSpec, theta_c: float, theta_amb: float, k_ins: float, k_jac: float,
                      Z=0.21, E=3.94, g=0.60, R=None):
    """IEC 60287 rating of a single cable in free air (no solar). Returns dict."""
    De = cable.outer_diameter
    h = Z / De**g + E
    T1, T3 = iec60287_T1_T3(cable, k_ins, k_jac)
    dtheta = theta_c - theta_amb
    KA = np.pi * De * h * (T1 + T3)
    s = 2.0
    for _ in range(200):
        s_new = (dtheta / (1 + KA * s)) ** 0.25
        if abs(s_new - s) < 1e-12:
            break
        s = s_new
    dtheta_s = s**4
    T4 = 1.0 / (np.pi * De * h * s)
    if R is None:
        R = float(cable.resistance_per_length(theta_c))
    I = np.sqrt(dtheta / (R * (T1 + T3 + T4)))
    return dict(I=I, T1=T1, T3=T3, T4=T4, h=h, dtheta_s=dtheta_s, R=R, W=I**2 * R)


def iec60287_conductor_temperature(cable: CableSpec, I: float, theta_amb: float, k_ins: float, k_jac: float, **kw):
    """Inverse of :func:`iec60287_free_air`: conductor temperature for a given current."""
    from scipy.optimize import brentq

    return brentq(lambda th: iec60287_free_air(cable, th, theta_amb, k_ins, k_jac, **kw)["I"] - I,
                  theta_amb + 1e-6, theta_amb + 400)


def iec60949_K(sigma_c: float, beta: float, rho20: float) -> float:
    """K in A s^0.5 / mm^2."""
    return np.sqrt(sigma_c * (beta + 20.0) / rho20) * 1e-6


def iec60949_final_temperature(I: float, t: float, S_mm2: float, theta_i: float, K: float = 226.0,
                               beta: float = 234.5) -> float:
    return (theta_i + beta) * np.exp(I**2 * t / (K**2 * S_mm2**2)) - beta


def iec60949_current(theta_i: float, theta_f: float, t: float, S_mm2: float, K: float = 226.0,
                     beta: float = 234.5) -> float:
    return K * S_mm2 * np.sqrt(np.log((theta_f + beta) / (theta_i + beta)) / t)
