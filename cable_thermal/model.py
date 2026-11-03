"""Radial finite-difference / TDMA model of an insulated and jacketed cable.

Faithful implementation of Section 2 of Riba & Llauradó, *Materials* 2022, 15,
6814 (Eqs. 2-27):

* the conductor (material A), inner insulation (B) and outer jacket (C) are
  discretized along the radius with N_A, N_B, N_C elements of thickness
  dx_A = R_A/N_A, dx_B = t_B/N_B, dx_C = t_C/N_C (Figs. 4-7);
* node 1 is the conductor centre, node N_A+1 the conductor-insulation
  boundary, node N_A+N_B+1 the insulation-jacket boundary and node
  N_A+N_B+N_C+1 the jacket-air boundary (total N_A+N_B+N_C+1 nodes);
* each node owns the half elements on both sides: areas S_i,A, S_i,B, S_i,C
  and face perimeters p_{i-1}, p_{i+1} exactly as Eqs. (6), (9), (12), (15),
  (18), (21), (24);
* backward (implicit) Euler in time; the Joule source I_i^2 r_i(T_i^j), the
  properties and the radiation term are evaluated at the previous time level j,
  the convection term h p (T^{j+1} - T_air) is implicit (Eqs. 5-27);
* the tri-diagonal system a_i T_{i-1} + b_i T_i + c_i T_{i+1} = d_i (Eqs. 3-4)
  is solved with the Thomas algorithm (TDMA).

Things the paper leaves open and how they are resolved (all configurable):

* temperature at which k is evaluated on a face -> mean of the two nodes
  (keeps a_{i+1} = c_i, i.e. the scheme is conservative);
* J_i = I_i/S_i with which total area -> the conductor region (radius R_A) is
  homogenized so that the total Joule loss is exactly I^2 r(T) with
  r = rho_e/S_eff (Eq. 2 / Table 1); ``joule_area='geometric'`` uses
  rho_e/(pi R_A^2) instead;
* conductor heat capacity -> homogenized density m/(pi R_A^2) so that
  sum(rho_A S_i,A c_p) = m c_p (consistent with Eq. 1 and Table 1 mass m).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence, Tuple, Union

import numpy as np

from .cables import CableSpec
from .heat_transfer import KELVIN, SIGMA, h_natural

Number = Union[int, float]


# ---------------------------------------------------------------------------
# Thomas algorithm (TDMA)
# ---------------------------------------------------------------------------
def tdma(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Solve a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i  (a_0 = c_{n-1} = 0)."""
    n = len(d)
    cp = np.empty(n)
    dp = np.empty(n)
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m if i < n - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = np.empty(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


# ---------------------------------------------------------------------------
# Load / ambient profiles
# ---------------------------------------------------------------------------
class StepProfile:
    """Piecewise-constant profile from ((value, duration), ...) (paper Tables 2 and 4)."""

    def __init__(self, steps: Sequence[Tuple[float, float]], after: Optional[float] = None):
        self.values = np.array([s[0] for s in steps], dtype=float)
        self.durations = np.array([s[1] for s in steps], dtype=float)
        self.edges = np.concatenate([[0.0], np.cumsum(self.durations)])
        self.after = self.values[-1] if after is None else after

    @property
    def total_duration(self) -> float:
        return float(self.edges[-1])

    def __call__(self, t):
        t = np.asarray(t, dtype=float)
        idx = np.searchsorted(self.edges, t, side="right") - 1
        vals = np.where(idx >= len(self.values), self.after, self.values[np.clip(idx, 0, len(self.values) - 1)])
        return vals if vals.ndim else float(vals)

    @property
    def breakpoints(self):
        return self.edges


def as_function(x) -> Callable[[float], float]:
    if callable(x):
        return x
    if isinstance(x, (list, tuple)) and len(x) and isinstance(x[0], (list, tuple)):
        return StepProfile(x)
    v = float(x)
    return lambda t: v


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
@dataclass
class ModelOptions:
    nodes_per_material: Tuple[int, int, int] = (9, 9, 9)   # paper: nine nodes per material
    dt: float = 10.0                                       # paper: 10 s
    convection_model: str = "ieee738"   # 'ieee738' (corrected Eq. 25) or 'paper' (as printed)
    elevation: float = 0.0              # H in Eq. (26) [m] (not given in the paper)
    joule_area: str = "effective"       # 'effective' -> r = rho/S_eff ; 'geometric' -> rho/(pi R_A^2)
    radiation: str = "explicit"         # 'explicit' (paper Eq. 23) or 'linearized' (Newton, implicit)
    face_k: str = "mean"                # k on a face evaluated at mean of the 2 nodes ('mean') or upwind node ('left')
    convection_off: bool = False        # verification switches
    radiation_off: bool = False
    adiabatic: bool = False             # no heat leaves the outer surface
    joule_scale: float = 1.0            # sensitivity/calibration multiplier on I^2 r (1.0 = paper)
    h_scale: float = 1.0                # sensitivity/calibration multiplier on h (1.0 = paper)
    rad_scale: float = 1.0              # sensitivity/calibration multiplier on radiation (1.0 = paper)
    capacity_scale: float = 1.0         # sensitivity/calibration multiplier on rho*c_p of all layers


@dataclass
class SimulationResult:
    t: np.ndarray                 # [s] (n_t,)
    T: np.ndarray                 # [°C] (n_t, n_nodes)
    x: np.ndarray                 # node radii [m]
    node_material: list
    probes: dict                  # index of the paper's measuring points
    current: np.ndarray           # [A] applied in step ending at t
    ambient: np.ndarray           # [°C]
    joule: np.ndarray             # [W/m]
    p_conv: np.ndarray            # [W/m]
    p_rad: np.ndarray             # [W/m]
    stored: np.ndarray            # [W/m] rate of heat storage (sum C dT/dt)
    cable: CableSpec = field(repr=False)
    options: ModelOptions = field(repr=False)

    def probe(self, name: str) -> np.ndarray:
        return self.T[:, self.probes[name]]

    @property
    def T_conductor_insulation(self):
        return self.probe("conductor_insulation")

    @property
    def T_insulation_jacket(self):
        return self.probe("insulation_jacket")

    @property
    def T_jacket_air(self):
        return self.probe("jacket_air")

    @property
    def T_center(self):
        return self.probe("center")

    def energy_balance_error(self) -> float:
        """max |stored - (Joule - Pc - Pr)| [W/m] - must be ~1e-12 (discrete conservation)."""
        return float(np.max(np.abs(self.stored[1:] - (self.joule[1:] - self.p_conv[1:] - self.p_rad[1:]))))

    def at(self, time: float) -> np.ndarray:
        j = int(np.argmin(np.abs(self.t - time)))
        return self.T[j]

    def to_csv(self, path: str):
        cols = ["t_s", "I_A", "T_air_C", "T_center_C", "T_cond_ins_C", "T_ins_jacket_C", "T_jacket_air_C",
                "joule_W_per_m", "P_conv_W_per_m", "P_rad_W_per_m"]
        data = np.column_stack([self.t, self.current, self.ambient, self.T_center, self.T_conductor_insulation,
                                self.T_insulation_jacket, self.T_jacket_air, self.joule, self.p_conv, self.p_rad])
        np.savetxt(path, data, delimiter=",", header=",".join(cols), comments="", fmt="%.6g")


class CableThermalModel:
    """Paper's FDM/TDMA cable model (see module docstring)."""

    def __init__(self, cable: CableSpec, options: Optional[ModelOptions] = None, **kw):
        self.cable = cable
        self.opt = options or ModelOptions()
        for k, v in kw.items():
            setattr(self.opt, k, v)
        self._build_mesh()

    # ------------------------------------------------------------------ mesh
    def _build_mesh(self):
        c = self.cable
        NA, NB, NC = self.opt.nodes_per_material
        dxA, dxB, dxC = c.R_A / NA, c.insulation_thickness / NB, c.jacket_thickness / NC
        x = np.concatenate([
            np.arange(NA + 1) * dxA,                     # nodes 1..NA+1   (centre .. A/B boundary)
            c.R_A + np.arange(1, NB + 1) * dxB,          # nodes NA+2..NA+NB+1
            c.R_B + np.arange(1, NC + 1) * dxC,          # nodes ..NA+NB+NC+1
        ])
        x[NA] = c.R_A
        x[NA + NB] = c.R_B
        x[-1] = c.R_C
        n = len(x)
        # segment s joins node s and s+1 ; material 0=A, 1=B, 2=C
        seg_mat = np.array([0] * NA + [1] * NB + [2] * NC)
        seg_dx = np.array([dxA] * NA + [dxB] * NB + [dxC] * NC)
        # face perimeter p (Eqs. 6, 9, 12, 15, 18, 21): p = 2 pi (x_s + 0.5 dx)
        self.face_p = 2.0 * np.pi * (x[:-1] + 0.5 * seg_dx)
        self.seg_mat, self.seg_dx = seg_mat, seg_dx
        # half-element areas owned by each node, per material (S_i,A, S_i,B, S_i,C)
        S = np.zeros((n, 3))
        for s in range(n - 1):
            xm = x[s] + 0.5 * seg_dx[s]
            S[s, seg_mat[s]] += np.pi * (xm**2 - x[s] ** 2)          # right half of node s
            S[s + 1, seg_mat[s]] += np.pi * (x[s + 1] ** 2 - xm**2)  # left half of node s+1
        self.S = S
        self.x = x
        self.n = n
        self.dx = (dxA, dxB, dxC)
        self.node_material = ["A"] * NA + ["A|B"] + ["B"] * (NB - 1) + ["B|C"] + ["C"] * (NC - 1) + ["C|air"]
        self.probes = {"center": 0, "conductor_insulation": NA, "insulation_jacket": NA + NB, "jacket_air": n - 1}
        # Joule share of each node: S_i,A / (pi R_A^2)
        self.joule_share = S[:, 0] / c.conductor_geometric_area
        self.rho = np.array([c.conductor_density_effective, c.insulation_material.density, c.jacket_material.density])
        self.mats = [c.conductor_material, c.insulation_material, c.jacket_material]
        if self.opt.joule_area == "effective":
            self._r_scale = 1.0
        elif self.opt.joule_area == "geometric":
            self._r_scale = c.effective_cross_section / c.conductor_geometric_area
        else:
            raise ValueError("joule_area must be 'effective' or 'geometric'")

    # ------------------------------------------------------------ properties
    def heat_capacity(self, T):
        """C_i = sum_m rho_m S_i,m c_p,m(T_i)  [J/(m K)] (left-hand sides of Eqs. 5-23)."""
        C = np.zeros_like(T)
        for m in range(3):
            C += self.rho[m] * self.S[:, m] * self.mats[m].cp(T)
        return self.opt.capacity_scale * C

    def face_conductance(self, T):
        """G_s = k_m(T_face) p_s / dx_m  [W/(m K)] (the k p/dx terms of Eqs. 5-27)."""
        if self.opt.face_k == "mean":
            Tf = 0.5 * (T[:-1] + T[1:])
        else:
            Tf = T[:-1]
        k = np.empty(self.n - 1)
        for m in range(3):
            sel = self.seg_mat == m
            k[sel] = self.mats[m].k(Tf[sel])
        return k * self.face_p / self.seg_dx

    def joule_sources(self, T, I):
        """I_i^2 r_i(T_i) for every node [W/m]; sums to I^2 r(T) for a uniform T."""
        r = self.cable.resistance_per_length(T) * self._r_scale
        return self.opt.joule_scale * I**2 * r * self.joule_share

    def surface_h(self, Ts, Tair):
        if self.opt.convection_off or self.opt.adiabatic:
            return 0.0
        return self.opt.h_scale * float(h_natural(Ts, Tair, self.cable.outer_diameter, self.opt.elevation,
                                                  self.opt.convection_model))

    def surface_radiation(self, Ts, Tair):
        if self.opt.radiation_off or self.opt.adiabatic:
            return 0.0
        p = 2.0 * np.pi * self.cable.R_C
        return self.opt.rad_scale * self.cable.emissivity * SIGMA * p * ((Ts + KELVIN) ** 4 - (Tair + KELVIN) ** 4)

    def radiation_slope(self, Ts):
        if self.opt.radiation_off or self.opt.adiabatic:
            return 0.0
        p = 2.0 * np.pi * self.cable.R_C
        return self.opt.rad_scale * 4.0 * self.cable.emissivity * SIGMA * p * (Ts + KELVIN) ** 3

    # ---------------------------------------------------------- assemble
    def _coefficients(self, T, I, Tair, dt):
        """TDMA coefficients a, b, c, d of Eqs. (7), (10), (13), (16), (19), (22), (27)."""
        G = self.face_conductance(T)
        C = self.heat_capacity(T)
        cap = C / dt if np.isfinite(dt) else np.zeros_like(C)
        a = np.zeros(self.n)
        c = np.zeros(self.n)
        a[1:] = -G
        c[:-1] = -G
        b = cap.copy()
        b[1:] += G
        b[:-1] += G
        q = self.joule_sources(T, I)
        d = cap * T + q
        # last node (Eq. 27)
        Ts = T[-1]
        h = self.surface_h(Ts, Tair)
        p = 2.0 * np.pi * self.cable.R_C
        b[-1] += h * p
        d[-1] += h * p * Tair
        prad = self.surface_radiation(Ts, Tair)
        if self.opt.radiation == "linearized" or not np.isfinite(dt):
            s = self.radiation_slope(Ts)
            b[-1] += s
            d[-1] += -prad + s * Ts
        else:
            d[-1] += -prad
        return a, b, c, d, C, q, h, p, prad

    # ---------------------------------------------------------- transient
    def simulate(self, current, ambient=None, t_end: Optional[float] = None, T_init=None,
                 dt: Optional[float] = None, record_every: int = 1) -> SimulationResult:
        """Transient / dynamic simulation.

        current : A, callable I(t), or ((I, duration), ...) step table
        ambient : °C, callable T_air(t), or step table (default: cable.ambient_temperature)
        """
        dt = float(dt or self.opt.dt)
        Ifun = as_function(current)
        if ambient is None:
            ambient = self.cable.ambient_temperature
        Tafun = as_function(ambient)
        if t_end is None:
            if isinstance(Ifun, StepProfile):
                t_end = Ifun.total_duration
            else:
                raise ValueError("t_end required")
        nsteps = int(round(t_end / dt))
        if T_init is None:
            T = np.full(self.n, float(Tafun(0.0)))
        elif np.isscalar(T_init):
            T = np.full(self.n, float(T_init))
        else:
            T = np.array(T_init, dtype=float)
        rec_t, rec_T, rec_I, rec_Ta, rec_q, rec_pc, rec_pr, rec_st = ([0.0], [T.copy()], [float(Ifun(0.0))],
                                                                     [float(Tafun(0.0))], [np.nan], [np.nan],
                                                                     [np.nan], [np.nan])
        t = 0.0
        for j in range(nsteps):
            tm = t + 0.5 * dt
            I = float(Ifun(tm))
            Tair = float(Tafun(t + dt))
            a, b, c, d, C, q, h, p, prad = self._coefficients(T, I, Tair, dt)
            Tn = tdma(a, b, c, d)
            if self.opt.radiation == "linearized":
                prad_used = prad + self.radiation_slope(T[-1]) * (Tn[-1] - T[-1])
            else:
                prad_used = prad
            stored = float(np.sum(C * (Tn - T)) / dt)
            t = (j + 1) * dt
            if (j + 1) % record_every == 0 or j == nsteps - 1:
                rec_t.append(t)
                rec_T.append(Tn.copy())
                rec_I.append(I)
                rec_Ta.append(Tair)
                rec_q.append(float(q.sum()))
                rec_pc.append(h * p * (Tn[-1] - Tair))
                rec_pr.append(prad_used)
                rec_st.append(stored)
            T = Tn
        return SimulationResult(np.array(rec_t), np.array(rec_T), self.x.copy(), self.node_material,
                                dict(self.probes), np.array(rec_I), np.array(rec_Ta), np.array(rec_q),
                                np.array(rec_pc), np.array(rec_pr), np.array(rec_st), self.cable, self.opt)

    # ---------------------------------------------------------- steady state
    def steady_state(self, current: float, ambient: Optional[float] = None, tol: float = 1e-9,
                     max_iter: int = 500, T_guess=None) -> np.ndarray:
        """Steady state (dt -> inf) of the same discrete equations, solved by Picard
        iteration with Newton-linearized radiation. Returns node temperatures [°C]."""
        Tair = self.cable.ambient_temperature if ambient is None else float(ambient)
        T = np.full(self.n, Tair + 10.0) if T_guess is None else np.array(T_guess, dtype=float)
        for it in range(max_iter):
            a, b, c, d, *_rest = self._coefficients(T, current, Tair, np.inf)
            Tn = tdma(a, b, c, d)
            if np.max(np.abs(Tn - T)) < tol:
                return Tn
            T = 0.5 * T + 0.5 * Tn if it > 50 else Tn
        raise RuntimeError("steady state did not converge")

    def steady_heat_terms(self, current: float, ambient: Optional[float] = None) -> dict:
        Tair = self.cable.ambient_temperature if ambient is None else float(ambient)
        T = self.steady_state(current, Tair)
        Ts = T[-1]
        h = self.surface_h(Ts, Tair)
        p = 2 * np.pi * self.cable.R_C
        return {"T": T, "joule": float(self.joule_sources(T, current).sum()), "p_conv": h * p * (Ts - Tair),
                "p_rad": float(self.surface_radiation(Ts, Tair)), "h_conv": h}

    def ampacity(self, T_max: Optional[float] = None, ambient: Optional[float] = None, location: str = "center") -> float:
        """Current giving T_max at the chosen node in steady state (current-temperature relationship).

        The upper bracket is grown gradually: far above the rating the linear r(T) law can
        lead to thermal runaway (no steady state), which is treated as 'too hot'."""
        from scipy.optimize import brentq

        T_max = self.cable.max_continuous_temperature if T_max is None else T_max
        idx = self.probes[location]

        def f(I):
            try:
                return self.steady_state(I, ambient)[idx] - T_max
            except (RuntimeError, FloatingPointError):
                return np.inf

        lo, hi = 0.0, 50.0
        with np.errstate(all="ignore"):
            while f(hi) < 0:
                lo, hi = hi, hi * 1.5
                if hi > 1e6:
                    raise RuntimeError("no ampacity bracket found")
            return brentq(f, lo, hi, xtol=1e-6)


def simulate_bare_conductor_eq1(cable: CableSpec, current, ambient=None, t_end=None, dt=1.0,
                                convection_model="ieee738", elevation=0.0):
    """Paper Eq. (1) (CIGRE/IEEE 738 lumped heat balance of a BARE conductor) with
    Eq. (2) resistance, solved with explicit Euler. Used only for comparison."""
    Ifun = as_function(current)
    Tafun = as_function(cable.ambient_temperature if ambient is None else ambient)
    if t_end is None:
        t_end = Ifun.total_duration
    D = cable.conductor_diameter
    m = cable.conductor_mass
    cp = cable.conductor_material.cp
    n = int(round(t_end / dt))
    T = np.empty(n + 1)
    T[0] = Tafun(0.0)
    for j in range(n):
        tj = j * dt
        Ta = Tafun(tj)
        Pj = Ifun(tj + 0.5 * dt) ** 2 * cable.resistance_per_length(T[j])
        h = h_natural(T[j], Ta, D, elevation, convection_model)
        Pc = h * np.pi * D * (T[j] - Ta)
        Pr = cable.emissivity * SIGMA * np.pi * D * ((T[j] + KELVIN) ** 4 - (Ta + KELVIN) ** 4)
        T[j + 1] = T[j] + dt * (Pj - Pc - Pr) / (m * cp(T[j]))
    return np.arange(n + 1) * dt, T
