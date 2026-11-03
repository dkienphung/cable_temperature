"""Independent high-resolution reference solver (verification only).

Same physics as the paper (radial conduction with temperature-dependent
properties, Joule source with r(T), natural convection + radiation at the
surface) but a different discretization, written independently of
:mod:`cable_thermal.model`:

* finite-volume method of lines on a fine grid;
* conductance between neighbouring nodes inside a layer from the exact
  steady cylindrical solution, G = 2 pi k / ln(r_{i+1}/r_i) (G = 2 pi k r_{1/2}/dr
  for the first segment that touches r = 0);
* radiation and convection fully implicit;
* time integration with SciPy's stiff Radau IIA (5th order) with tight tolerances.

Used to measure the spatial/temporal discretization error of the paper's scheme
(9 nodes per material, dt = 10 s).
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags

from .cables import CableSpec
from .heat_transfer import KELVIN, SIGMA, h_natural
from .model import StepProfile, as_function


class ReferenceSolver:
    def __init__(self, cable: CableSpec, n=(60, 60, 60), convection_model="ieee738", elevation=0.0):
        self.cable = cable
        self.conv = convection_model
        self.elev = elevation
        c = cable
        NA, NB, NC = n
        x = np.concatenate([np.linspace(0, c.R_A, NA + 1), np.linspace(c.R_A, c.R_B, NB + 1)[1:],
                            np.linspace(c.R_B, c.R_C, NC + 1)[1:]])
        self.x = x
        self.n = len(x)
        self.seg_mat = np.array([0] * NA + [1] * NB + [2] * NC)
        xm = 0.5 * (x[:-1] + x[1:])
        # control-volume areas by material
        self.V = np.zeros((self.n, 3))
        for s in range(self.n - 1):
            self.V[s, self.seg_mat[s]] += np.pi * (xm[s] ** 2 - x[s] ** 2)
            self.V[s + 1, self.seg_mat[s]] += np.pi * (x[s + 1] ** 2 - xm[s] ** 2)
        self.geo = np.empty(self.n - 1)
        for s in range(self.n - 1):
            if x[s] == 0.0:
                self.geo[s] = 2 * np.pi * xm[s] / (x[s + 1] - x[s])
            else:
                self.geo[s] = 2 * np.pi / np.log(x[s + 1] / x[s])
        self.rho = np.array([c.conductor_density_effective, c.insulation_material.density, c.jacket_material.density])
        self.mats = [c.conductor_material, c.insulation_material, c.jacket_material]
        self.share = self.V[:, 0] / c.conductor_geometric_area
        self.p = 2 * np.pi * c.R_C
        self.probes = {"center": 0, "conductor_insulation": NA, "insulation_jacket": NA + NB, "jacket_air": self.n - 1}

    def rhs(self, t, T, Ifun, Tafun):
        c = self.cable
        Tf = 0.5 * (T[:-1] + T[1:])
        k = np.empty(self.n - 1)
        for m in range(3):
            sel = self.seg_mat == m
            k[sel] = self.mats[m].k(Tf[sel])
        G = k * self.geo
        flux = G * (T[:-1] - T[1:])           # heat flow node s -> s+1 [W/m]
        net = np.zeros(self.n)
        net[:-1] -= flux
        net[1:] += flux
        I = Ifun(t)
        net += I**2 * c.resistance_per_length(T) * self.share
        Ta = Tafun(t)
        Ts = T[-1]
        h = float(h_natural(Ts, Ta, c.outer_diameter, self.elev, self.conv))
        net[-1] -= h * self.p * (Ts - Ta) + c.emissivity * SIGMA * self.p * ((Ts + KELVIN) ** 4 - (Ta + KELVIN) ** 4)
        C = np.zeros(self.n)
        for m in range(3):
            C += self.rho[m] * self.V[:, m] * self.mats[m].cp(T)
        return net / C

    def simulate(self, current, ambient=None, t_end=None, t_eval=None, rtol=1e-8, atol=1e-8):
        Ifun = as_function(current)
        Tafun = as_function(self.cable.ambient_temperature if ambient is None else ambient)
        if t_end is None:
            t_end = Ifun.total_duration
        bps = list(Ifun.breakpoints) if isinstance(Ifun, StepProfile) else [0.0, t_end]
        bps = sorted(set([b for b in bps if b < t_end] + [0.0, t_end]))
        if t_eval is None:
            t_eval = np.linspace(0, t_end, int(t_end // 10) + 1)
        T = np.full(self.n, float(Tafun(0.0)))
        jac_sp = diags([np.ones(self.n - 1), np.ones(self.n), np.ones(self.n - 1)], [-1, 0, 1])
        ts, Ts = [], []
        for a, b in zip(bps[:-1], bps[1:]):
            # the step value on (a, b) is evaluated inside the interval
            sel = t_eval[(t_eval >= a) & (t_eval <= b)]
            sol = solve_ivp(self.rhs, (a, b), T, method="Radau", t_eval=sel, rtol=rtol, atol=atol,
                            jac_sparsity=jac_sp, args=(lambda t, a=a, b=b: Ifun(0.5 * (a + b)) if isinstance(Ifun, StepProfile) else Ifun(t), Tafun))
            if not sol.success:
                raise RuntimeError(sol.message)
            T = sol.y[:, -1]
            for tt, yy in zip(sol.t, sol.y.T):
                if ts and abs(tt - ts[-1]) < 1e-9:
                    continue
                ts.append(tt)
                Ts.append(yy)
        return np.array(ts), np.array(Ts)


def analytic_steady_state_constant_props(cable: CableSpec, q: float, T_surface: float, kA, kB, kC):
    """Exact steady temperatures for constant k and uniform heat generation q [W/m]
    in the conductor: returns (T_center, T_A|B, T_B|C)."""
    T_BC = T_surface + q / (2 * np.pi * kC) * np.log(cable.R_C / cable.R_B)
    T_AB = T_BC + q / (2 * np.pi * kB) * np.log(cable.R_B / cable.R_A)
    T_0 = T_AB + q / (4 * np.pi * kA)
    return T_0, T_AB, T_BC
