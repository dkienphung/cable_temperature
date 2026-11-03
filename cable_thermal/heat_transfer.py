"""Surface heat-transfer terms of the paper (Eqs. 23-27) and their check.

Natural convection - paper Eq. (25), taken from IEEE Std 738-2012 [26]::

    h = 3.645/pi * rho_air^0.25 * D^-0.25 * (T - T_air)^0.25          (as printed)

IEEE Std 738-2012, Eq. (5c), gives the natural-convection heat loss per unit length::

    q_cn = 3.645 * rho_f^0.5 * D0^0.75 * (Ts - Ta)^1.25                 [W/m]

Dividing by the surface area per unit length (pi*D) and by (Ts - Ta) gives::

    h = q_cn / (pi D (Ts - Ta)) = 3.645/pi * rho_f^0.5 * D^-0.25 * (Ts - Ta)^0.25

so the exponent of the air density printed in Eq. (25) (0.25) is wrong; it must be
0.5. Both versions are implemented (``convection_model='paper'`` reproduces the
printed formula, ``'ieee738'`` the corrected one).

Air density - paper Eq. (26), IEEE Std 738-2012 Eq. (14a)::

    rho_air = (1.293 - 1.525e-4 H + 6.379e-9 H^2) / (1 + 0.00367 T_film)

The paper states the temperatures in kelvin, but the formula (0.00367 = 1/272.5)
requires T_film = (T_air + T) / 2 in degrees Celsius. With kelvin the density
would be ~0.62 kg/m³ at room temperature instead of ~1.2 kg/m³. The code always
uses °C (``film_temperature_in_kelvin=True`` exists only to quantify the error).

Radiation - paper Eq. (23)/(27)::

    P_r = eps * sigma * pi * D * (T_s^4 - T_air^4)  with temperatures in K, eps = 0.85
"""
from __future__ import annotations

import numpy as np

SIGMA = 5.670374419e-8  # Stefan–Boltzmann constant [W/(m² K⁴)] (CODATA 2018)
KELVIN = 273.15


def air_density_ieee738(T_surface_C, T_air_C, elevation_m: float = 0.0, film_temperature_in_kelvin: bool = False):
    """Eq. (26) / IEEE 738-2012 Eq. (14a). Temperatures in °C. Returns kg/m³."""
    T_film = 0.5 * (np.asarray(T_surface_C) + np.asarray(T_air_C))
    if film_temperature_in_kelvin:  # literal reading of the paper's text (wrong)
        T_film = T_film + KELVIN
    H = elevation_m
    return (1.293 - 1.525e-4 * H + 6.379e-9 * H**2) / (1.0 + 0.00367 * T_film)


def h_natural(T_surface_C, T_air_C, D: float, elevation_m: float = 0.0, model: str = "ieee738"):
    """Natural-convection heat-transfer coefficient h [W/(m² K)] for a horizontal
    cylinder of outer diameter D [m].

    model = 'paper'   -> Eq. (25) exactly as printed (rho^0.25)
    model = 'ieee738' -> corrected Eq. (25) consistent with IEEE 738-2012 (rho^0.5)
    """
    dT = np.abs(np.asarray(T_surface_C, dtype=float) - T_air_C)
    rho = air_density_ieee738(T_surface_C, T_air_C, elevation_m)
    if model == "paper":
        exp_rho = 0.25
    elif model == "ieee738":
        exp_rho = 0.5
    else:
        raise ValueError(f"unknown convection model {model!r}")
    return 3.645 / np.pi * rho**exp_rho * D**-0.25 * dT**0.25


def q_convection_ieee738(T_surface_C, T_air_C, D: float, elevation_m: float = 0.0):
    """IEEE 738-2012 Eq. (5c), natural convection loss per unit length [W/m] (reference)."""
    rho = air_density_ieee738(T_surface_C, T_air_C, elevation_m)
    dT = np.asarray(T_surface_C, dtype=float) - T_air_C
    return 3.645 * rho**0.5 * D**0.75 * np.sign(dT) * np.abs(dT) ** 1.25


def q_radiation(T_surface_C, T_air_C, D: float, emissivity: float = 0.85):
    """Radiated power per unit length [W/m], Eq. (23): eps*sigma*pi*D*(Ts^4 - Ta^4)."""
    Ts = np.asarray(T_surface_C, dtype=float) + KELVIN
    Ta = np.asarray(T_air_C, dtype=float) + KELVIN
    return emissivity * SIGMA * np.pi * D * (Ts**4 - Ta**4)
