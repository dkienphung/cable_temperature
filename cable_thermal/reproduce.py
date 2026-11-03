"""Reproduce the results/figures of Riba & Llauradó (2022) and the extra studies.

Run:  python -m cable_thermal reproduce --out results
"""
from __future__ import annotations

import json
import os
from dataclasses import asdict
from importlib import resources

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from .cables import CABLE_70, CABLE_150, CableSpec
from .heat_transfer import air_density_ieee738, h_natural, q_radiation
from .iec import iec60287_free_air, iec60949_current, iec60949_final_temperature, iec60949_K
from .model import CableThermalModel, ModelOptions, simulate_bare_conductor_eq1
from .properties import digitized, material_at_reference, pvc, xlpe
from .reference import ReferenceSolver

# palette (reference categorical palette, slots chosen to mirror the paper's colour code)
ORANGE, BLUE, GREY, VIOLET, AQUA, RED, YELLOW = "#eb6834", "#2a78d6", "#6b6a65", "#6250d6", "#1baf7a", "#e34948", "#eda100"
INK, INK2 = "#0b0b0b", "#52514e"
PROBES = [("conductor_insulation", "Conductor-insulation", ORANGE),
          ("insulation_jacket", "Insulation-jacket", BLUE),
          ("jacket_air", "Jacket-air", GREY)]

REF = json.loads(resources.files("cable_thermal.data").joinpath("paper_figures_digitized.json").read_text())


def _style():
    plt.rcParams.update({
        "figure.dpi": 100, "savefig.dpi": 200, "font.size": 10, "axes.titlesize": 11,
        "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2,
        "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
        "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.8,
        "legend.frameon": False, "lines.linewidth": 1.8, "axes.axisbelow": True,
    })


def _paper_curve(fig_key, probe):
    d = REF[fig_key]
    key = {"jacket_air": "jacket_air_measured"}.get(probe, probe)
    return np.array(d[key])


def _rms_bias(result, fig_key):
    out = {}
    allerr = []
    for probe, _, _ in PROBES:
        d = _paper_curve(fig_key, probe)
        e = np.interp(d[:, 0], result.t, result.probe(probe)) - d[:, 1]
        out[probe] = {"rms_K": float(np.sqrt(np.mean(e**2))), "mean_K": float(e.mean()), "max_abs_K": float(np.abs(e).max()),
                      "mean_rel_pct": float(np.mean(100 * e / d[:, 1])), "end_model_C": float(result.probe(probe)[-1]),
                      "end_paper_C": float(d[-1, 1])}
        allerr.append(e)
    e = np.concatenate(allerr)
    out["all"] = {"rms_K": float(np.sqrt(np.mean(e**2))), "mean_K": float(e.mean())}
    return out


def _step_labels(ax, cable, result, y_offset=0.6):
    edges = np.concatenate([[0], np.cumsum([s[1] for s in cable.test_steps])])
    for (I, _), a, b in zip(cable.test_steps, edges[:-1], edges[1:]):
        tm = 0.5 * (a + b)
        y = np.interp(tm, result.t, result.T_conductor_insulation)
        ax.text(tm, y + y_offset, f"{I:g} A", ha="center", va="bottom", fontsize=9, fontweight="bold", color=INK)


def fig_temperature_rise(result, cable, fig_key, path, title, extra=None, legend_loc="lower right"):
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.2, 7.6), height_ratios=[1.5, 1], sharex=True)
    for probe, label, col in PROBES:
        ax.plot(result.t, result.probe(probe), color=col, label=f"{label} (this model)")
        d = _paper_curve(fig_key, probe)
        src = "paper, measured" if probe == "jacket_air" else "paper Fig., meas.+sim."
        ax.plot(d[:, 0], d[:, 1], color=col, ls=(0, (4, 2)), lw=1.2, alpha=0.9, label=f"{label} ({src}, digitized)")
    if extra:
        for t, T, lab, col, ls in extra:
            ax.plot(t, T, color=col, ls=ls, lw=1.4, label=lab)
    _step_labels(ax, cable, result)
    ax.set_ylabel("Temperature [°C]")
    ax.set_title(title, loc="left")
    ax.legend(fontsize=7.5, loc=legend_loc, ncol=1)
    for probe, label, col in PROBES:
        d = _paper_curve(fig_key, probe)
        Tm = np.interp(d[:, 0], result.t, result.probe(probe))
        ax2.plot(d[:, 0], 100 * (d[:, 1] - Tm) / d[:, 1], color=col, label=label, lw=1.4)
    ax2.axhline(0, color=INK2, lw=0.8)
    ax2.set_ylabel("Relative difference [%]\n(paper − model)/paper")
    ax2.set_xlabel("Time [s]")
    ax2.legend(fontsize=8, ncol=3, loc="lower left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_properties(path_cp, path_k):
    T = np.linspace(15, 75, 241)
    for path, keys, ylab, title in [
        (path_cp, [("copper_cp", "Copper"), ("xlpe_cp", "XLPE"), ("pvc_cp", "PVC")], "c$_p$ [J/(kg K)]",
         "Specific heat capacity (digitized from paper Fig. 2)"),
        (path_k, [("copper_k", "Copper"), ("xlpe_k", "XLPE"), ("pvc_k", "PVC")], "k [W/(m K)]",
         "Thermal conductivity (digitized from paper Fig. 3)"),
    ]:
        fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
        for ax, (key, name), col in zip(axs, keys, [ORANGE, BLUE, GREY]):
            f = digitized(key)
            Tr = T[(T >= f.T[0])]
            ax.plot(Tr, f(Tr), color=col)
            ax.plot(f.T, f.values, "o", ms=3, color=col, alpha=0.6)
            ax.set_title(name, loc="left")
            ax.set_xlabel("Temperature [°C]")
            ax.set_ylabel(ylab)
        fig.suptitle(title, x=0.01, ha="left", fontsize=11)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)


def fig_discretization(model, path):
    fig, ax = plt.subplots(figsize=(10, 1.9))
    c = model.cable
    for (a, b, col, name) in [(0, c.R_A, "#f6c26b", "Conductor (A)"), (c.R_A, c.R_B, "#8fbbe8", "Insulation XLPE (B)"),
                              (c.R_B, c.R_C, "#a9a8a3", "Jacket PVC (C)")]:
        ax.axvspan(a * 1e3, b * 1e3, color=col, alpha=0.8, lw=0)
        ax.text(0.5 * (a + b) * 1e3, -0.55, name, ha="center", va="top", fontsize=9)
    ax.plot(model.x * 1e3, np.zeros_like(model.x), "o", color=INK, ms=5)
    for xf in (model.x[:-1] + 0.5 * model.seg_dx):
        ax.axvline(xf * 1e3, color=INK2, ls=":", lw=0.8)
    for i in [0, model.probes["conductor_insulation"], model.probes["insulation_jacket"], model.n - 1]:
        ax.annotate(f"i={i + 1}", (model.x[i] * 1e3, 0), xytext=(0, 9), textcoords="offset points", ha="center", fontsize=8)
    ax.set_ylim(-1.2, 1)
    ax.set_yticks([])
    ax.grid(False)
    ax.set_xlabel("Radius [mm]   (dots: nodes, dotted lines: element faces x_i ± Δx/2)")
    ax.set_title(f"Radial discretization of the {c.name} cable, N_A = N_B = N_C = {model.opt.nodes_per_material[0]} "
                 f"({model.n} nodes) - paper Figs. 4-7", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_radial_profiles(result, model, path, times):
    fig, ax = plt.subplots(figsize=(7, 5))
    NA, NB, NC = model.opt.nodes_per_material
    nodes = np.arange(1, model.n + 1)
    ax.axvspan(1, NA + 1, color="#f6c26b", alpha=0.35, lw=0)
    ax.axvspan(NA + 1, NA + NB + 1, color="#8fbbe8", alpha=0.35, lw=0)
    ax.axvspan(NA + NB + 1, model.n, color="#a9a8a3", alpha=0.35, lw=0)
    for lab, xm in [("CONDUCTOR", (1 + NA + 1) / 2), ("INSULATION", (2 * NA + NB + 2) / 2), ("JACKET", (NA + NB + 1 + model.n) / 2)]:
        ax.text(xm, 1.01, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=9, color=INK2)
    for t in times:
        T = result.at(t)
        ax.plot(nodes, T, color=INK, lw=1.2)
        ax.text(nodes[NA + NB // 2], T[NA + NB // 2] + 0.25, f"t = {t:g} s", fontsize=8, ha="center", color=INK)
    ax.set_xlabel("Node position")
    ax.set_ylabel("Temperature [°C]")
    ax.set_xlim(0, model.n + 1)
    ax.set_title("Temporal evolution of the radial temperature (cf. paper Fig. 10a)", loc="left", fontsize=10, pad=18)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_heat_balance(result, path):
    fig, ax = plt.subplots(figsize=(7, 5))
    t = result.t[1:]
    ax.plot(t, result.joule[1:], color=INK, lw=2.2, label="Joule heating I²r (model)")
    ax.plot(t, result.p_conv[1:] + result.p_rad[1:], color=INK, ls=(0, (5, 3)), lw=2, label="P$_c$ + P$_r$ (model)")
    ax.plot(t, result.p_conv[1:], color=ORANGE, label="Convection heat loss P$_c$ (model)")
    ax.plot(t, result.p_rad[1:], color=BLUE, label="Radiation heat loss P$_r$ (model)")
    d = REF["fig10b_cable70_heat_terms"]
    pc, pr = np.array(d["Pc"]), np.array(d["Pr"])
    ax.plot(pc[:, 0], pc[:, 1], color=ORANGE, ls=":", lw=1.6, label="P$_c$ (paper Fig. 10b, digitized)")
    ax.plot(pr[:, 0], pr[:, 1], color=BLUE, ls=":", lw=1.6, label="P$_r$ (paper Fig. 10b, digitized)")
    lv = d["joule_levels_read_by_eye"]
    edges = np.concatenate([[0], np.cumsum([s[1] for s in CABLE_70.test_steps])])
    for k, a, b in zip(["120A", "145A", "170A", "185A", "212A"], edges[:-1], edges[1:]):
        ax.plot([a, b], lv[k], color=GREY, ls=":", lw=2.2, label="I²r (paper Fig. 10b, read)" if k == "120A" else None)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Per unit length heat terms [W/m]")
    ax.set_title("Components of the heat balance, 70 mm² (cf. paper Fig. 10b)", loc="left", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def calibrate(cable, fig_key, names=("capacity_scale", "h_scale")):
    """Diagnostic inverse analysis: multipliers that make the model match the paper curves."""
    def resid(p):
        r = CableThermalModel(cable, **dict(zip(names, p))).simulate(cable.test_steps)
        out = []
        for probe, _, _ in PROBES:
            d = _paper_curve(fig_key, probe)
            out.append(np.interp(d[:, 0], r.t, r.probe(probe)) - d[:, 1])
        return np.concatenate(out)
    sol = least_squares(resid, np.ones(len(names)), bounds=(0.2, 5.0))
    return dict(zip(names, map(float, sol.x))), float(np.sqrt(np.mean(sol.fun**2)))


def run_all(out="results", quick=False):
    _style()
    figdir = os.path.join(out, "figures")
    datadir = os.path.join(out, "data")
    os.makedirs(figdir, exist_ok=True)
    os.makedirs(datadir, exist_ok=True)
    F = lambda n: os.path.join(figdir, n)
    summary = {}

    # ---------------------------------------------------------------- Figs. 2-3, 4
    fig_properties(F("fig02_specific_heat.png"), F("fig03_thermal_conductivity.png"))
    m70 = CableThermalModel(CABLE_70)
    fig_discretization(m70, F("fig04_discretization.png"))
    summary["cables"] = {"cable70": CABLE_70.summary(), "cable150": CABLE_150.summary()}

    # ---------------------------------------------------------------- Fig. 9 / 10 / 11 (as published)
    runs = {}
    for cable, key in [(CABLE_70, "fig9_cable70"), (CABLE_150, "fig11_cable150")]:
        for conv in ("paper", "ieee738"):
            r = CableThermalModel(cable, convection_model=conv).simulate(cable.test_steps)
            runs[(cable.name, conv)] = r
            summary.setdefault("comparison_with_paper", {})[f"{cable.name}/{conv}"] = _rms_bias(r, key)
        runs[(cable.name, "ieee738")].to_csv(os.path.join(datadir, f"{cable.name}_test_sequence.csv"))

    r70 = runs[("cable70", "ieee738")]
    r150 = runs[("cable150", "ieee738")]
    fig_temperature_rise(r70, CABLE_70, "fig9_cable70", F("fig09_cable70_temperature_rise.png"),
                         "Conductor of 70 mm²: model (Table 1 data, corrected Eq. 25) vs paper Fig. 9")
    fig_temperature_rise(r150, CABLE_150, "fig11_cable150", F("fig11_cable150_temperature_rise.png"),
                         "Conductor of 150 mm²: model (Table 3 data, corrected Eq. 25) vs paper Fig. 11",
                         legend_loc="upper left")
    fig_radial_profiles(r70, m70, F("fig10a_cable70_radial_profiles.png"), times=list(range(0, 8001, 1000)))
    fig_heat_balance(r70, F("fig10b_cable70_heat_balance.png"))
    summary["heat_balance_cable70_end"] = {"joule": float(r70.joule[-1]), "p_conv": float(r70.p_conv[-1]),
                                           "p_rad": float(r70.p_rad[-1]), "T_surface": float(r70.T_jacket_air[-1])}
    summary["energy_balance_error_W_per_m"] = {k[0] + "/" + k[1]: r.energy_balance_error() for k, r in runs.items()}

    # ---------------------------------------------------------------- Eq. (25) printed vs corrected + diagnostic fit
    cal = {}
    # the 150 mm2 test (3100 s) is far from steady state, so h is not identifiable there -> fit capacity only
    for cable, key, names in [(CABLE_70, "fig9_cable70", ("capacity_scale", "h_scale")),
                              (CABLE_150, "fig11_cable150", ("capacity_scale",))]:
        p, rms = calibrate(cable, key, names)
        cal[cable.name] = {"factors": p, "rms_K": rms}
    summary["diagnostic_calibration"] = cal
    rcal = CableThermalModel(CABLE_70, **cal["cable70"]["factors"]).simulate(CABLE_70.test_steps)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    d = _paper_curve("fig9_cable70", "conductor_insulation")
    ax.plot(d[:, 0], d[:, 1], color=INK, ls=(0, (4, 2)), lw=1.4, label="Paper Fig. 9 (digitized)")
    ax.plot(runs[("cable70", "paper")].t, runs[("cable70", "paper")].T_conductor_insulation, color=RED,
            label="Model, Eq. (25) as printed (ρ$^{0.25}$)")
    ax.plot(r70.t, r70.T_conductor_insulation, color=ORANGE, label="Model, Eq. (25) corrected (ρ$^{0.5}$, IEEE 738)")
    f = cal["cable70"]["factors"]
    ax.plot(rcal.t, rcal.T_conductor_insulation, color=VIOLET, ls="-.",
            label=f"Diagnostic fit: heat capacity ×{f['capacity_scale']:.2f}, h ×{f['h_scale']:.2f}")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Conductor-insulation temperature [°C]")
    ax.set_title("70 mm²: effect of the Eq. (25) correction and diagnostic inverse fit", loc="left", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(F("fig09b_cable70_formula_variants.png"))
    plt.close(fig)

    # ---------------------------------------------------------------- verification vs reference solver + convergence
    ver = {}
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for cable, col in [(CABLE_70, ORANGE), (CABLE_150, BLUE)]:
        ref = ReferenceSolver(cable, n=(40, 40, 40))
        tr, Tr = ref.simulate(cable.test_steps)
        r = runs[(cable.name, "ieee738")]
        diff = np.interp(tr, r.t, r.T_conductor_insulation) - Tr[:, ref.probes["conductor_insulation"]]
        axs[0].plot(tr, diff, color=col, lw=1.3, label=f"{cable.name}: paper scheme (9 nodes, Δt=10 s) − reference")
        conv = []
        for N, dt in [(3, 40), (5, 20), (9, 10), (18, 5), (36, 2.5)]:
            rr = CableThermalModel(cable, nodes_per_material=(N, N, N), dt=dt).simulate(cable.test_steps)
            e = np.max(np.abs(np.interp(tr, rr.t, rr.T_conductor_insulation) - Tr[:, ref.probes["conductor_insulation"]]))
            conv.append((N, dt, float(e)))
        ver[cable.name] = {"max_abs_diff_paper_scheme_K": float(np.abs(diff).max()), "convergence": conv}
        axs[1].loglog([c[0] for c in conv], [c[2] for c in conv], "o-", color=col, label=cable.name)
    axs[0].set_xlabel("Time [s]")
    axs[0].set_ylabel("ΔT conductor-insulation [K]")
    axs[0].set_title("Paper FDM/TDMA vs independent Radau reference (121 nodes)", loc="left", fontsize=10)
    axs[0].legend(fontsize=8)
    axs[1].set_xlabel("Nodes per material N (Δt refined with N: 40, 20, 10, 5, 2.5 s)")
    axs[1].set_ylabel("max |T − T$_{ref}$| [K]")
    axs[1].set_title("Grid/time-step convergence", loc="left", fontsize=10)
    axs[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(F("fig12_verification_reference_solver.png"))
    plt.close(fig)
    summary["verification"] = ver

    # ---------------------------------------------------------------- steady-state current-temperature relationship
    st = {}
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
    for ax, cable in zip(axs, [CABLE_70, CABLE_150]):
        m = CableThermalModel(cable)
        Ta = cable.ambient_temperature
        Imax = m.ampacity(90.0)
        I = np.linspace(0, 1.15 * Imax, 40)
        T = np.array([m.steady_state(i) for i in I])
        ax.plot(I, T[:, 0], color=ORANGE, label="Conductor (model, steady state)")
        ax.plot(I, T[:, -1], color=GREY, label="Jacket surface (model)")
        # IEC 60287 with the paper's material data (k at 75 °C, the end of the paper's curves)
        kB, kC = float(cable.insulation_material.k(75)), float(cable.jacket_material.k(75))
        Ti = np.linspace(Ta + 1, 100, 40)
        Iiec = [iec60287_free_air(cable, th, Ta, kB, kC)["I"] for th in Ti]
        ax.plot(Iiec, Ti, color=VIOLET, ls="--", label="IEC 60287 (k from paper curves)")
        Iiec_std = [iec60287_free_air(cable, th, Ta, 1 / 3.5, 1 / 5.0)["I"] for th in Ti]
        ax.plot(Iiec_std, Ti, color=AQUA, ls=":", label="IEC 60287 (IEC ρ$_T$: XLPE 3.5, PVC 5.0)")
        ax.axhline(90, color=RED, lw=0.9)
        ax.annotate(f"{Imax:.0f} A @ 90 °C", (Imax, 90), xytext=(-70, 10), textcoords="offset points", fontsize=8,
                    arrowprops=dict(arrowstyle="->", color=INK2))
        # paper test currents
        for (Ik, _) in cable.test_steps:
            ax.axvline(Ik, color="#d8d7d2", lw=0.8, zorder=0)
        iec90 = iec60287_free_air(cable, 90.0, Ta, kB, kC)
        iec90s = iec60287_free_air(cable, 90.0, Ta, 1 / 3.5, 1 / 5.0)
        st[cable.name] = {"ambient_C": Ta, "ampacity_90C_model_A": Imax, "ampacity_90C_IEC60287_paper_k_A": iec90["I"],
                          "ampacity_90C_IEC60287_IEC_rhoT_A": iec90s["I"],
                          "note": "property curves are clamped above 75 degC"}
        ax.set_xlabel("Current [A rms]")
        ax.set_ylabel("Temperature [°C]")
        ax.set_title(f"{cable.name}: steady-state current-temperature (T$_{{air}}$ = {Ta:g} °C)", loc="left", fontsize=10)
        ax.legend(fontsize=7.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(F("fig13_steady_state_current_temperature.png"))
    plt.close(fig)
    summary["steady_state"] = st

    # ---------------------------------------------------------------- short circuit
    sc = {}
    cable = CABLE_70
    m = CableThermalModel(cable)
    Iload = st["cable70"]["ampacity_90C_model_A"]
    T0 = m.steady_state(Iload)
    alpha = cable.alpha
    beta_p = 1.0 / alpha - cable.T0
    cp_ref = float(cable.conductor_material.cp(75.0))
    K_p = np.sqrt(cable.conductor_mass * cp_ref / (cable.resistance_per_length(cable.T0) * alpha)) / (cable.effective_cross_section * 1e6)
    t_sc = 1.0
    Isc = iec60949_current(T0[0], 250.0, t_sc, 70.0, K_p, beta_p)
    prof = lambda t: Isc if t < t_sc else Iload
    rsc = m.simulate(prof, t_end=60.0, dt=0.005, T_init=T0, record_every=4)
    r10 = m.simulate(prof, t_end=60.0, dt=10.0, T_init=T0)
    rsc_long = m.simulate(prof, t_end=3600.0, dt=0.5, T_init=T0)
    tt = np.linspace(0, t_sc, 101)
    iec_curve = iec60949_final_temperature(Isc, tt, 70.0, T0[0], K_p, beta_p)
    # IEC 60949 non-adiabatic factor eps = 1 + X sqrt(t/S) + Y t/S (Cu conductor, XLPE insulation: X=0.41, Y=0.12)
    eps_na = 1 + 0.41 * np.sqrt(t_sc / 70.0) + 0.12 * t_sc / 70.0
    T_nonadiabatic = iec60949_final_temperature(Isc / eps_na, t_sc, 70.0, T0[0], K_p, beta_p)
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.0))
    ax = axs[0]
    for probe, label, col in PROBES:
        ax.plot(rsc.t, rsc.probe(probe), color=col, label=f"{label} (model, Δt = 5 ms)")
    ax.plot(rsc.t, rsc.T_center, color=INK, lw=1, label="Conductor centre (model)")
    ax.plot(tt, iec_curve, color=VIOLET, ls="--", label="IEC 60949 adiabatic (same Cu data)")
    ax.plot([t_sc], [T_nonadiabatic], "D", color=AQUA, ms=7, label=f"IEC 60949 non-adiabatic (ε = {eps_na:.3f})")
    ax.plot(r10.t, r10.T_conductor_insulation, "s", color=RED, ms=5, label="Paper time step Δt = 10 s")
    ax.set_xlim(0, 20)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [°C]")
    ax.set_title(f"70 mm²: {Isc / 1e3:.2f} kA, {t_sc:g} s short circuit from {T0[0]:.0f} °C", loc="left", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    ax = axs[1]
    for probe, label, col in PROBES:
        ax.plot(rsc_long.t / 60, rsc_long.probe(probe), color=col, label=label)
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Temperature [°C]")
    ax.set_title(f"Recovery at load current {Iload:.0f} A", loc="left", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(F("fig14_short_circuit.png"))
    plt.close(fig)
    sc["cable70"] = {"I_load_A": Iload, "T_initial_C": float(T0[0]), "I_sc_A": float(Isc), "t_sc_s": t_sc,
                     "K_paper_data": float(K_p), "beta_paper_data": beta_p,
                     "K_IEC_copper": float(iec60949_K(3.45e6, 234.5, 1.7241e-8)),
                     "T_max_model_center_C": float(rsc.T_center.max()),
                     "T_max_model_conductor_insulation_C": float(rsc.T_conductor_insulation.max()),
                     "T_final_IEC60949_adiabatic_C": float(iec_curve[-1]),
                     "IEC60949_non_adiabatic_factor": float(eps_na),
                     "T_final_IEC60949_non_adiabatic_C": float(T_nonadiabatic),
                     "T_max_paper_dt10_C": float(r10.T_conductor_insulation.max())}
    summary["short_circuit"] = sc

    # ---------------------------------------------------------------- dynamic: variable current and ambient
    cable = CABLE_70
    m = CableThermalModel(cable)
    hours = 48
    t_end = hours * 3600
    load = lambda t: 120 + 80 * np.clip(np.sin(2 * np.pi * (t / 86400 - 0.3)), 0, None) + (60 if 30 * 3600 <= t < 32 * 3600 else 0)
    amb = lambda t: 19 + 6 * np.sin(2 * np.pi * (t / 86400 - 0.375))
    rd = m.simulate(load, ambient=amb, t_end=t_end, dt=30.0, record_every=4)
    rd_const = m.simulate(load, ambient=19.0, t_end=t_end, dt=30.0, record_every=4)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6.2), sharex=True, height_ratios=[1, 2])
    ax1.plot(rd.t / 3600, rd.current, color=INK, lw=1.4)
    ax1.set_ylabel("Current [A]")
    ax1b = ax1.twinx()
    ax1b.plot(rd.t / 3600, rd.ambient, color=AQUA, lw=1.4)
    ax1b.set_ylabel("T$_{air}$ [°C]", color=AQUA)
    ax1b.spines["right"].set_visible(True)
    ax1.set_title("Dynamic problem, 70 mm²: daily load cycle + 2 h overload, daily ambient cycle", loc="left", fontsize=10)
    for probe, label, col in PROBES:
        ax2.plot(rd.t / 3600, rd.probe(probe), color=col, label=f"{label} (variable ambient)")
    ax2.plot(rd_const.t / 3600, rd_const.T_conductor_insulation, color=ORANGE, ls=":", label="Conductor-insulation (T$_{air}$ = 19 °C)")
    ax2.plot(rd.t / 3600, rd.ambient, color=AQUA, lw=1, label="Ambient")
    ax2.set_xlabel("Time [h]")
    ax2.set_ylabel("Temperature [°C]")
    ax2.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(F("fig15_dynamic_current_ambient.png"))
    plt.close(fig)
    summary["dynamic"] = {"T_max_conductor_C": float(rd.T_conductor_insulation.max()),
                          "T_max_conductor_const_ambient_C": float(rd_const.T_conductor_insulation.max())}

    # ---------------------------------------------------------------- sensitivity (end of test, 70 mm²)
    base = r70.T_conductor_insulation[-1]
    variants = [
        ("Eq. (25) as printed (ρ^0.25)", CABLE_70, dict(convection_model="paper")),
        ("α = 0.00393 1/K (Cu @ 20 °C)", CABLE_70.with_(alpha=0.00393), {}),
        ("ρ_e = 1.7241e-8 Ω·m (IEC 60028)", CABLE_70.with_(resistivity_T0=1.7241e-8), {}),
        ("J over π R_A² (joule_area=geometric)", CABLE_70, dict(joule_area="geometric")),
        ("ε = 0.95", CABLE_70.with_(emissivity=0.95), {}),
        ("ε = 0.75", CABLE_70.with_(emissivity=0.75), {}),
        ("Elevation H = 280 m (Terrassa)", CABLE_70, dict(elevation=280.0)),
        ("XLPE & PVC density +20 %", CABLE_70.with_(insulation_material=xlpe(920 * 1.2), jacket_material=pvc(1400 * 1.2)), {}),
        ("XLPE & PVC density −20 %", CABLE_70.with_(insulation_material=xlpe(920 * 0.8), jacket_material=pvc(1400 * 0.8)), {}),
        ("Conductor mass = 8960·S_eff", CABLE_70.with_(conductor_mass_per_length=8960 * 70e-6), {}),
        ("Constant properties @20 °C", CABLE_70.with_(conductor_material=material_at_reference(CABLE_70.conductor_material),
                                                      insulation_material=material_at_reference(CABLE_70.insulation_material),
                                                      jacket_material=material_at_reference(CABLE_70.jacket_material)), {}),
        ("Radiation linearized (implicit)", CABLE_70, dict(radiation="linearized")),
        ("Fig. 10b I²r level (joule ×0.90)", CABLE_70, dict(joule_scale=0.90)),
    ]
    sens = []
    for lab, c, kw in variants:
        r = CableThermalModel(c, **kw).simulate(c.test_steps)
        sens.append((lab, float(r.T_conductor_insulation[-1] - base), float(_rms_bias(r, "fig9_cable70")["all"]["rms_K"])))
    summary["sensitivity_cable70_end_of_test"] = {"base_C": float(base), "variants": sens}
    fig, ax = plt.subplots(figsize=(8, 5))
    labs = [s[0] for s in sens][::-1]
    vals = [s[1] for s in sens][::-1]
    ax.barh(labs, vals, color=[ORANGE if v > 0 else BLUE for v in vals], height=0.6)
    for y, v in enumerate(vals):
        ax.text(v + (0.05 if v >= 0 else -0.05), y, f"{v:+.2f} K", va="center", ha="left" if v >= 0 else "right", fontsize=8)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_xlabel(f"Change of conductor-insulation temperature at t = 8600 s [K]  (base {base:.2f} °C)")
    ax.set_title("70 mm²: sensitivity to formula choices and uncertain data", loc="left", fontsize=10)
    lim = max(abs(min(vals)), abs(max(vals))) * 1.35
    ax.set_xlim(-lim, lim)
    fig.tight_layout()
    fig.savefig(F("fig16_sensitivity.png"))
    plt.close(fig)

    # ---------------------------------------------------------------- Eq. (1) bare conductor vs insulated cable
    tb, Tb = simulate_bare_conductor_eq1(CABLE_70, CABLE_70.test_steps, dt=1.0)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(r70.t, r70.T_center, color=ORANGE, label="Insulated + jacketed cable (FDM model), conductor")
    ax.plot(tb, Tb, color=INK, ls="--", label="Same conductor bare, lumped Eq. (1) (CIGRE / IEEE 738)")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Conductor temperature [°C]")
    ax.set_title("70 mm²: insulated cable vs bare-conductor heat balance Eq. (1)", loc="left", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(F("fig17_bare_conductor_eq1.png"))
    plt.close(fig)
    summary["bare_conductor_eq1_end_C"] = float(Tb[-1])

    # ---------------------------------------------------------------- h vs ΔT (Eq. 25 printed vs corrected) and Eq. 26 units
    fig, ax = plt.subplots(figsize=(7, 4))
    dT = np.linspace(0.5, 60, 200)
    D = CABLE_70.outer_diameter
    Ta = 19.0
    ax.plot(dT, h_natural(Ta + dT, Ta, D, model="ieee738"), color=ORANGE, label="h, Eq. (25) corrected (ρ$^{0.5}$)")
    ax.plot(dT, h_natural(Ta + dT, Ta, D, model="paper"), color=RED, ls="--", label="h, Eq. (25) as printed (ρ$^{0.25}$)")
    hr = q_radiation(Ta + dT, Ta, D, 0.85) / (np.pi * D * dT)
    ax.plot(dT, hr, color=BLUE, label="Equivalent radiation coefficient (ε = 0.85)")
    rhoK = air_density_ieee738(Ta + dT, Ta, film_temperature_in_kelvin=True)
    ax.plot(dT, 3.645 / np.pi * rhoK**0.5 * D**-0.25 * dT**0.25, color=GREY, ls=":",
            label="h with Eq. (26) evaluated in K (text of paper) - wrong")
    ax.set_xlabel("Surface temperature rise T$_s$ − T$_{air}$ [K]")
    ax.set_ylabel("Heat-transfer coefficient [W/(m² K)]")
    ax.set_title(f"Surface heat-transfer coefficients, D = {D * 1e3:.1f} mm", loc="left", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(F("fig18_heat_transfer_coefficients.png"))
    plt.close(fig)

    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    write_results_md(summary, os.path.join(out, "RESULTS.md"))
    return summary


def write_results_md(s, path):
    c70, c150 = s["comparison_with_paper"]["cable70/ieee738"], s["comparison_with_paper"]["cable150/ieee738"]
    p70, p150 = s["comparison_with_paper"]["cable70/paper"], s["comparison_with_paper"]["cable150/paper"]
    L = []
    w = L.append
    w("# Results - cable current-temperature model (Riba & Llauradó, Materials 2022, 15, 6814)\n")
    w("All numbers below are produced by `python -m cable_thermal reproduce` (file generated automatically).\n")
    w("Paper curves were digitized from the published figures (±0.15 °C reading accuracy). For the "
      "conductor-insulation and insulation-jacket traces the paper's measured and simulated curves overlap and "
      "the digitized trace is their median; the jacket-air trace is the measured one.\n")
    w("## 1. Test sequences of the paper (Figs. 9 and 11)\n")
    w("Model = paper Eqs. (2)-(27) with the cable data of Tables 1-4, 9 nodes per material, Δt = 10 s.\n")
    for name, c, p in [("70 mm² (Table 1/2, Fig. 9)", c70, p70), ("150 mm² (Table 3/4, Fig. 11)", c150, p150)]:
        w(f"### {name}\n")
        w("| Measuring point | paper end [°C] | model end, Eq. 25 corrected [°C] | model end, Eq. 25 as printed [°C] | mean error [K] | RMS error [K] | mean rel. error [%] |")
        w("|---|---|---|---|---|---|---|")
        for probe, label, _ in PROBES:
            q, qp = c[probe], p[probe]
            w(f"| {label} | {q['end_paper_C']:.2f} | {q['end_model_C']:.2f} | {qp['end_model_C']:.2f} | {q['mean_K']:+.2f} | {q['rms_K']:.2f} | {q['mean_rel_pct']:+.1f} |")
        w("")
    w("Discrete energy balance error of the TDMA solution (max over all steps, W/m): "
      + ", ".join(f"{k}: {v:.1e}" for k, v in s["energy_balance_error_W_per_m"].items()) + "\n")
    hb = s["heat_balance_cable70_end"]
    w(f"Heat balance at t = 8600 s (70 mm², model): I²r = {hb['joule']:.2f} W/m, P_c = {hb['p_conv']:.2f} W/m, "
      f"P_r = {hb['p_rad']:.2f} W/m, surface {hb['T_surface']:.2f} °C. Paper Fig. 10b: I²r ≈ 10.8, P_c ≈ 6.3, P_r ≈ 4.4 W/m.\n")
    w("## 2. Diagnostic inverse analysis\n")
    w("Least-squares multipliers on the total heat capacity and on the convection coefficient that make the model "
      "reproduce the digitized paper curves (diagnosis only - not used elsewhere):\n")
    w("| Cable | heat-capacity factor | h factor | RMS after fit [K] |")
    w("|---|---|---|---|")
    for k, v in s["diagnostic_calibration"].items():
        hs = f"{v['factors']['h_scale']:.2f}" if "h_scale" in v["factors"] else "1 (not fitted)"
        w(f"| {k} | {v['factors']['capacity_scale']:.2f} | {hs} | {v['rms_K']:.2f} |")
    w("")
    w("## 3. Numerical verification\n")
    for k, v in s["verification"].items():
        conv = ", ".join(f"N={n} (Δt={dt:g} s): {e:.4f} K" for n, dt, e in v["convergence"])
        w(f"* {k}: paper scheme vs independent 121-node Radau reference, max |ΔT| = {v['max_abs_diff_paper_scheme_K']:.4f} K. Convergence: {conv}.")
    w("")
    w("## 4. Steady-state current-temperature relationship (ampacity at 90 °C)\n")
    w("| Cable | T_air [°C] | model [A] | IEC 60287, k from paper curves [A] | IEC 60287, IEC ρ_T [A] |")
    w("|---|---|---|---|---|")
    for k, v in s["steady_state"].items():
        w(f"| {k} | {v['ambient_C']:g} | {v['ampacity_90C_model_A']:.0f} | {v['ampacity_90C_IEC60287_paper_k_A']:.0f} | {v['ampacity_90C_IEC60287_IEC_rhoT_A']:.0f} |")
    w("\nThe paper's property curves end at 75 °C and are held constant above (documented limitation).\n")
    scv = s["short_circuit"]["cable70"]
    w("## 5. Short circuit (70 mm²)\n")
    w(f"Load {scv['I_load_A']:.0f} A (steady {scv['T_initial_C']:.1f} °C), then {scv['I_sc_A'] / 1e3:.2f} kA for {scv['t_sc_s']:g} s "
      f"(the current that IEC 60949 predicts gives 250 °C adiabatically with the same copper data, K = {scv['K_paper_data']:.1f} A·s^0.5/mm², "
      f"β = {scv['beta_paper_data']:.1f} K; standard IEC copper K = {scv['K_IEC_copper']:.1f}).\n")
    w(f"* model peak conductor temperature: {scv['T_max_model_center_C']:.1f} °C (centre), {scv['T_max_model_conductor_insulation_C']:.1f} °C (conductor-insulation boundary)")
    w(f"* IEC 60949 adiabatic: {scv['T_final_IEC60949_adiabatic_C']:.1f} °C (upper bound: no heat leaves the conductor)")
    w(f"* IEC 60949 non-adiabatic (ε = 1 + X√(t/S) + Y t/S, Cu/XLPE X = 0.41, Y = 0.12 → ε = {scv['IEC60949_non_adiabatic_factor']:.3f}): "
      f"{scv['T_final_IEC60949_non_adiabatic_C']:.1f} °C - the model's heat flow into the XLPE gives practically the same result")
    w(f"* with the paper's Δt = 10 s the 1 s pulse is not resolved (peak {scv['T_max_paper_dt10_C']:.1f} °C) - short circuits need Δt ≪ t_sc.\n")
    w("## 6. Sensitivity of the 70 mm² end-of-test temperature\n")
    sb = s["sensitivity_cable70_end_of_test"]
    w(f"Base (Table 1 data, Eq. 25 corrected): {sb['base_C']:.2f} °C at t = 8600 s; paper: {c70['conductor_insulation']['end_paper_C']:.2f} °C.\n")
    w("| Variant | ΔT [K] | RMS vs paper [K] |")
    w("|---|---|---|")
    for lab, dT, rms in sb["variants"]:
        w(f"| {lab} | {dT:+.2f} | {rms:.2f} |")
    w("")
    w("## 7. Dynamic case and Eq. (1)\n")
    w(f"Peak conductor-insulation temperature over 48 h: {s['dynamic']['T_max_conductor_C']:.1f} °C "
      f"(variable ambient) vs {s['dynamic']['T_max_conductor_const_ambient_C']:.1f} °C (constant 19 °C).\n")
    w(f"Bare conductor with lumped Eq. (1): {s['bare_conductor_eq1_end_C']:.2f} °C at 8600 s.\n")
    with open(path, "w") as fh:
        fh.write("\n".join(L))
