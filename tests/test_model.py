"""Verification tests of the formulas and of the numerical implementation."""
import numpy as np
import pytest

from cable_thermal import CABLE_70, CABLE_150, CableThermalModel, tdma
from cable_thermal.heat_transfer import (air_density_ieee738, h_natural, q_convection_ieee738, q_radiation)
from cable_thermal.iec import (iec60287_free_air, iec60949_K, iec60949_final_temperature, iec60287_conductor_temperature)
from cable_thermal.properties import constant_material, digitized
from cable_thermal.reference import ReferenceSolver, analytic_steady_state_constant_props


def const_cable(cable, kA=398.0, kB=0.33, kC=0.14):
    return cable.with_(conductor_material=constant_material("cu", 8960, 390.0, kA),
                       insulation_material=constant_material("ins", 920, 2000.0, kB),
                       jacket_material=constant_material("jac", 1400, 1000.0, kC))


# ---------------------------------------------------------------- TDMA (Eqs. 3-4)
def test_tdma_matches_dense_solver():
    rng = np.random.default_rng(0)
    n = 28
    a, c = -rng.random(n), -rng.random(n)
    a[0] = c[-1] = 0.0
    b = 3 + rng.random(n)
    d = rng.random(n)
    A = np.diag(b) + np.diag(a[1:], -1) + np.diag(c[:-1], 1)
    assert np.allclose(tdma(a, b, c, d), np.linalg.solve(A, d), atol=1e-13)


# ---------------------------------------------------------------- mesh (Eqs. 6, 9, 12, 15, 18, 21, 24)
@pytest.mark.parametrize("cable", [CABLE_70, CABLE_150])
def test_mesh_matches_paper(cable):
    m = CableThermalModel(cable)
    NA, NB, NC = m.opt.nodes_per_material
    assert m.n == NA + NB + NC + 1 == 28            # 28 nodal elements as in Fig. 10a
    assert np.isclose(m.x[NA], cable.R_A) and np.isclose(m.x[NA + NB], cable.R_B) and np.isclose(m.x[-1], cable.R_C)
    dxA, dxB, dxC = m.dx
    # node 1 (Eq. 6)
    assert np.isclose(m.S[0, 0], np.pi * (0.5 * dxA) ** 2)
    assert np.isclose(m.face_p[0], 2 * np.pi * 0.5 * dxA)
    # generic conductor node (Eq. 9)
    i = 3
    assert np.isclose(m.S[i, 0], np.pi * ((m.x[i] + 0.5 * dxA) ** 2 - (m.x[i] - 0.5 * dxA) ** 2))
    # boundary node A|B (Eq. 12)
    i = NA
    assert np.isclose(m.S[i, 0], np.pi * (m.x[i] ** 2 - (m.x[i] - 0.5 * dxA) ** 2))
    assert np.isclose(m.S[i, 1], np.pi * ((m.x[i] + 0.5 * dxB) ** 2 - m.x[i] ** 2))
    # boundary node B|C (Eq. 18) and last node (Eq. 24)
    i = NA + NB
    assert np.isclose(m.S[i, 2], np.pi * ((m.x[i] + 0.5 * dxC) ** 2 - m.x[i] ** 2))
    assert np.isclose(m.S[-1, 2], np.pi * (cable.R_C ** 2 - (cable.R_C - 0.5 * dxC) ** 2))
    # areas add up to the real layer areas
    assert np.isclose(m.S[:, 0].sum(), np.pi * cable.R_A ** 2)
    assert np.isclose(m.S[:, 1].sum(), np.pi * (cable.R_B ** 2 - cable.R_A ** 2))
    assert np.isclose(m.S[:, 2].sum(), np.pi * (cable.R_C ** 2 - cable.R_B ** 2))
    assert np.isclose(m.joule_share.sum(), 1.0)


def test_total_joule_equals_I2r_and_mass_is_table_value():
    m = CableThermalModel(CABLE_70)
    T = np.full(m.n, 40.0)
    assert np.isclose(m.joule_sources(T, 200.0).sum(), 200.0 ** 2 * CABLE_70.resistance_per_length(40.0))
    # Eq. (2) with Table 1: r20 = 1.85e-8 / 70e-6
    assert np.isclose(CABLE_70.resistance_per_length(20.0), 1.85e-8 / 70e-6)
    cp = CABLE_70.conductor_material.cp(30.0)
    assert np.isclose((m.rho[0] * m.S[:, 0]).sum(), 0.584)
    assert np.isclose(m.heat_capacity(np.full(m.n, 30.0))[: 9].sum() + m.rho[0] * m.S[9, 0] * cp,
                      0.584 * cp)


# ---------------------------------------------------------------- conservation
@pytest.mark.parametrize("radiation", ["explicit", "linearized"])
def test_discrete_energy_conservation(radiation):
    r = CableThermalModel(CABLE_70, radiation=radiation).simulate(CABLE_70.test_steps)
    assert r.energy_balance_error() < 1e-8


# ---------------------------------------------------------------- steady state vs analytic
@pytest.mark.parametrize("N,tol", [(9, 0.02), (40, 0.002)])
def test_steady_state_against_analytic(N, tol):
    cable = const_cable(CABLE_70)
    cable0 = cable.with_(alpha=0.0)
    m = CableThermalModel(cable0, nodes_per_material=(N, N, N))
    I = 200.0
    T = m.steady_state(I)
    q = I ** 2 * cable0.resistance_per_length(20.0)
    T0, TAB, TBC = analytic_steady_state_constant_props(cable0, q, T[-1], 398.0, 0.33, 0.14)
    dT = T0 - T[-1]
    assert abs(T[m.probes["insulation_jacket"]] - TBC) < tol * dT
    assert abs(T[m.probes["conductor_insulation"]] - TAB) < tol * dT
    assert abs(T[0] - T0) < tol * dT
    # surface balance: q = Pc + Pr
    h = m.surface_h(T[-1], cable0.ambient_temperature)
    loss = h * 2 * np.pi * cable0.R_C * (T[-1] - cable0.ambient_temperature) + m.surface_radiation(T[-1], cable0.ambient_temperature)
    assert np.isclose(loss, q, rtol=1e-8)


def test_long_transient_reaches_steady_state():
    m = CableThermalModel(CABLE_70)
    Tss = m.steady_state(150.0)
    r = m.simulate(150.0, t_end=40000.0, dt=20.0)
    assert np.max(np.abs(r.T[-1] - Tss)) < 1e-3


# ---------------------------------------------------------------- short circuit vs IEC 60949 (adiabatic)
def test_adiabatic_short_circuit_matches_iec60949():
    rho20, beta, sigma_c = 1.7241e-8, 234.5, 3.45e6
    S = 70e-6
    alpha = 1.0 / (beta + 20.0)
    cp = sigma_c / 8960.0
    cable = CABLE_70.with_(resistivity_T0=rho20, alpha=alpha, effective_cross_section=S,
                           conductor_mass_per_length=8960.0 * S,
                           conductor_material=constant_material("cu", 8960, cp, 398.0),
                           insulation_material=constant_material("ins", 1e-9, 1.0, 1e-12),
                           jacket_material=constant_material("jac", 1e-9, 1.0, 1e-12))
    K = iec60949_K(sigma_c, beta, rho20)
    assert abs(K - 226.0) < 0.5
    I, t = 10_000.0, 1.0
    m = CableThermalModel(cable, adiabatic=True, nodes_per_material=(9, 9, 9))
    r = m.simulate(I, ambient=90.0, t_end=t, dt=1e-4, T_init=90.0, record_every=1000)
    theta_f = iec60949_final_temperature(I, t, 70.0, 90.0, K, beta)
    assert abs(r.T_center[-1] - theta_f) < 0.05


# ---------------------------------------------------------------- surface heat-transfer formulas
def test_corrected_eq25_equals_ieee738_heat_loss():
    D, Ts, Ta = CABLE_70.outer_diameter, 45.0, 19.0
    h = h_natural(Ts, Ta, D, model="ieee738")
    assert np.isclose(h * np.pi * D * (Ts - Ta), q_convection_ieee738(Ts, Ta, D))
    # the printed Eq. (25) differs by rho^0.25
    rho = air_density_ieee738(Ts, Ta)
    assert np.isclose(h_natural(Ts, Ta, D, model="paper") * rho ** 0.25, h)


def test_air_density_requires_celsius():
    assert np.isclose(air_density_ieee738(20.0, 20.0), 1.293 / (1 + 0.00367 * 20.0))
    assert 1.15 < air_density_ieee738(30.0, 20.0) < 1.25            # physical value
    assert air_density_ieee738(30.0, 20.0, film_temperature_in_kelvin=True) < 0.7   # literal K reading: wrong


def test_radiation_formula():
    D = 0.0147
    assert np.isclose(q_radiation(35.0, 19.0, D, 0.85),
                      0.85 * 5.670374419e-8 * np.pi * D * ((35 + 273.15) ** 4 - (19 + 273.15) ** 4))


# ---------------------------------------------------------------- IEC 60287 iteration self-consistency
def test_iec60287_free_air_consistency():
    res = iec60287_free_air(CABLE_70, 90.0, 30.0, k_ins=1 / 3.5, k_jac=1 / 5.0)
    W = res["W"]
    # surface rise from T4 and heat flow must agree with the iterated value
    assert np.isclose(W * res["T4"], res["dtheta_s"], rtol=1e-8)
    assert np.isclose(W * (res["T1"] + res["T3"] + res["T4"]), 60.0, rtol=1e-8)
    th = iec60287_conductor_temperature(CABLE_70, res["I"], 30.0, 1 / 3.5, 1 / 5.0)
    assert abs(th - 90.0) < 1e-6


# ---------------------------------------------------------------- independent reference solver
@pytest.mark.parametrize("cable", [CABLE_70, CABLE_150])
def test_paper_scheme_agrees_with_reference_solver(cable):
    ref = ReferenceSolver(cable, n=(30, 30, 30))
    t, T = ref.simulate(cable.test_steps)
    r = CableThermalModel(cable).simulate(cable.test_steps)
    for name in ("conductor_insulation", "insulation_jacket", "jacket_air"):
        assert np.max(np.abs(np.interp(t, r.t, r.probe(name)) - T[:, ref.probes[name]])) < 0.05


# ---------------------------------------------------------------- digitized properties
def test_digitized_properties_are_plausible():
    assert 385 < digitized("copper_cp")(20) < 389
    assert 1700 < digitized("xlpe_cp")(20) < 1800 and 3100 < digitized("xlpe_cp")(75) < 3200
    assert 940 < digitized("pvc_cp")(20) < 960
    assert 398 < digitized("copper_k")(20) < 399.5
    assert 0.31 < digitized("xlpe_k")(60) < 0.32
    assert 0.11 < digitized("pvc_k")(75) < 0.12


# ---------------------------------------------------------------- site elevation and indoor air movement
from cable_thermal.heat_transfer import (air_conductivity_ieee738, air_viscosity_ieee738, h_surface,
                                         q_forced_ieee738, reynolds_ieee738, wind_direction_factor)
from cable_thermal.model import ModelOptions


def test_paper_cables_use_lab_elevation_and_override():
    assert CABLE_70.elevation == CABLE_150.elevation == 280.0
    assert CableThermalModel(CABLE_70).elevation == 280.0
    assert CableThermalModel(CABLE_70, ModelOptions(elevation=0.0)).elevation == 0.0
    # Eq. (26): density at 280 m is lower than at sea level by the IEEE 738 polynomial ratio
    ratio = air_density_ieee738(30, 20, 280.0) / air_density_ieee738(30, 20, 0.0)
    assert np.isclose(ratio, (1.293 - 1.525e-4 * 280 + 6.379e-9 * 280 ** 2) / 1.293)
    hot_280 = CableThermalModel(CABLE_70).steady_state(212.0)[0]
    hot_0 = CableThermalModel(CABLE_70, elevation=0.0).steady_state(212.0)[0]
    assert 0.0 < hot_280 - hot_0 < 0.5          # thinner air -> slightly hotter


def test_ieee738_forced_convection_formulas():
    assert np.isclose(wind_direction_factor(90.0), 1.0)
    assert np.isclose(wind_direction_factor(0.0), 0.388)
    assert np.isclose(air_viscosity_ieee738(20.0), 1.458e-6 * 293 ** 1.5 / 403.4)
    assert 0.0255 < air_conductivity_ieee738(20.0) < 0.0260
    D, Ts, Ta, V = 0.0147, 40.0, 19.0, 0.5
    Re = reynolds_ieee738(Ts, Ta, D, V)
    kf = air_conductivity_ieee738(0.5 * (Ts + Ta))
    qc1 = (1.01 + 1.35 * Re ** 0.52) * kf * (Ts - Ta)
    qc2 = 0.754 * Re ** 0.6 * kf * (Ts - Ta)
    assert np.isclose(q_forced_ieee738(Ts, Ta, D, V), max(qc1, qc2))
    # h_surface = largest of natural and forced convection (IEEE 738)
    h = h_surface(Ts, Ta, D, air_speed=V)
    assert np.isclose(h, max(h_natural(Ts, Ta, D), max(qc1, qc2) / (np.pi * D * (Ts - Ta))))


def test_zero_air_speed_is_the_paper_still_air_model():
    D = CABLE_70.outer_diameter
    for Ts in (19.5, 30.0, 60.0):
        assert h_surface(Ts, 19.0, D, 280.0, air_speed=0.0) == h_natural(Ts, 19.0, D, 280.0)
    r0 = CableThermalModel(CABLE_70).simulate(CABLE_70.test_steps)
    r1 = CableThermalModel(CABLE_70, air_speed=0.0).simulate(CABLE_70.test_steps)
    assert np.array_equal(r0.T, r1.T)


def test_air_movement_cools_and_is_conservative():
    ends = [CableThermalModel(CABLE_70, air_speed=v).simulate(CABLE_70.test_steps) for v in (0.0, 0.2, 1.0)]
    T_end = [r.T_conductor_insulation[-1] for r in ends]
    assert T_end[0] > T_end[1] > T_end[2]
    assert all(r.energy_balance_error() < 1e-8 for r in ends)


def test_reference_solver_with_air_speed():
    ref = ReferenceSolver(CABLE_70, n=(30, 30, 30), air_speed=0.2)
    t, T = ref.simulate(CABLE_70.test_steps)
    r = CableThermalModel(CABLE_70, air_speed=0.2).simulate(CABLE_70.test_steps)
    assert np.max(np.abs(np.interp(t, r.t, r.T_conductor_insulation) - T[:, ref.probes["conductor_insulation"]])) < 0.05
