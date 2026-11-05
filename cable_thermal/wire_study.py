"""Wire test campaign: run the cable model on the whole wire database and compare with
manufacturer current ratings.

Run:  python -m cable_thermal wire-study --out results/wires
"""
from __future__ import annotations

import csv
import json
import math
import os
from collections import Counter, defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .heat_transfer import h_surface, q_radiation
from .model import CableThermalModel
from .reproduce import AQUA, BLUE, GREY, INK, INK2, ORANGE, RED, VIOLET, YELLOW, _style
from .wires import (all_alpha_wires, awg_sort_key, find_part, from_habia, from_rs, from_supervutron,
                    material_db, reference_ratings)

GREEN, MAGENTA = "#008300", "#e87ba4"
GROUPS = {  # display group -> (material keys, colour)
    "PVC": (("PVC",), ORANGE),
    "XLPE (+PVC jacket)": (("XLPE",), YELLOW),
    "Silicone rubber": (("silicone_rubber",), BLUE),
    "PTFE": (("PTFE",), VIOLET),
    "ETFE": (("ETFE",), AQUA),
    "Fluoroplastic (ThermoThin)": (("fluoropolymer_unspecified",), MAGENTA),
    "PVDF": (("PVDF",), GREEN),
    "mPPE (EcoWire)": (("mPPE",), RED),
    "EPR": (("EPR",), GREY),
}
AMBIENT = 30.0
CONVECTION = "churchill_chu"


def group_of(w) -> str:
    key = w.layers[0][0]
    for g, (keys, _) in GROUPS.items():
        if key in keys:
            return g
    return "other"


def time_constant(model: CableThermalModel, current: float, ambient: float) -> float:
    """Time for the conductor to reach 63.2 % of its steady temperature rise (step from ambient)."""
    Tss = model.steady_state(current, ambient)
    rise = Tss[0] - ambient
    C = model.heat_capacity(np.full(model.n, ambient)).sum()
    Ts = Tss[-1]
    D = model.cable.outer_diameter
    G = (float(model.surface_h(Ts, ambient)) * math.pi * D + q_radiation(Ts, ambient, D, model.cable.emissivity)
         / max(Ts - ambient, 1e-9))
    tau_est = C / G
    dt = tau_est / 300.0
    r = model.simulate(current, ambient=ambient, t_end=6 * tau_est, dt=dt)
    y = r.T_center - ambient
    j = int(np.argmax(y >= 0.632 * rise))
    if j == 0:
        return float("nan")
    return float(np.interp(0.632 * rise, y[j - 1:j + 1], r.t[j - 1:j + 1]))


def steady_T(model, I, Ta):
    return float(model.steady_state(I, Ta)[0])


# ----------------------------------------------------------------------------- campaign
def run_campaign(out_csv: str):
    rows = []
    for w in all_alpha_wires():
        row = {"part": w.part_number, "family": (w.raw or {}).get("family"), "awg": w.awg, "insulation": w.insulation_label,
               "group": group_of(w) if w.simulatable else "not simulated", "T_max_C": w.temp_max_C,
               "D_cond_mm": w.conductor_diameter * 1e3, "OD_mm": w.outer_diameter * 1e3,
               "wall_mm": (w.outer_diameter - w.conductor_diameter) / 2 * 1e3, "metal_area_mm2": w.metal_area * 1e6,
               "R20_ohm_per_km": w.r20 * 1e3, "R20_source": w.r20_source, "voltage_V": (w.raw or {}).get("voltage_V"),
               "simulated": w.simulatable, "notes": "; ".join(w.notes)}
        if w.simulatable and w.temp_max_C:
            c = w.to_cable(AMBIENT)
            m = CableThermalModel(c, convection_model=CONVECTION)
            m738 = CableThermalModel(c, convection_model="ieee738")
            Tmax = min(float(w.temp_max_C), 300.0)
            row["I_rated_A"] = m.ampacity(Tmax, AMBIENT)
            row["I_rated_A_eq25"] = m738.ampacity(Tmax, AMBIENT)
            row["I_105C_A"] = m.ampacity(105.0, AMBIENT)
            row["I_rise40K_A"] = m.ampacity(AMBIENT + 40.0, AMBIENT)
            row["tau_s"] = time_constant(m, row["I_rise40K_A"], AMBIENT)
            Tss = m.steady_state(row["I_rated_A"], AMBIENT)
            row["dT_insulation_at_rated_K"] = float(Tss[0] - Tss[-1])
            row["T_surface_at_rated_C"] = float(Tss[-1])
            row["current_density_A_mm2"] = row["I_rated_A"] / row["metal_area_mm2"]
        rows.append(row)
    keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if k in rows[0] else 99)
    with open(out_csv, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    return rows


# ----------------------------------------------------------------------------- reference comparisons
def run_comparisons():
    refs = reference_ratings()
    res = {}
    # Habia E (PTFE, 200 °C), ratings at 40 °C ambient
    meta = refs["habia_E_ptfe"]
    hab = []
    for row in meta["rows"]:
        w = from_habia(row, meta)
        out = {"type": row["type"], "awg": row["awg"], "rating_A": row["current_rating_A"]}
        for conv in ("churchill_chu", "ieee738"):
            m = CableThermalModel(w.to_cable(40.0), convection_model=conv)
            out[f"T_at_rating_{conv}"] = steady_T(m, row["current_rating_A"], 40.0)
            out[f"I_200C_{conv}"] = m.ampacity(200.0, 40.0)
        hab.append(out)
    # ambient derating predicted by the model for a 200 °C conductor
    der = {}
    for t in ("E 3007", "E 2219", "E 1419"):
        row = next(r for r in meta["rows"] if r["type"] == t)
        m = CableThermalModel(from_habia(row, meta).to_cable(40.0), convection_model=CONVECTION)
        I40 = m.ampacity(200.0, 40.0)
        der[t] = {int(Ta): m.ampacity(200.0, float(Ta)) / I40 for Ta in meta["ambient_correction"]}
    res["habia"] = {"rows": hab, "derating_model": der, "derating_published": meta["ambient_correction"]}
    # RS PRO PTFE Type A: rating = 40 K rise, ambient not stated -> 20 °C used (rise is almost independent)
    meta = refs["rs_ptfe_type_a"]
    rs = []
    for row in meta["rows"]:
        w = from_rs(row, meta)
        out = {"awg": row["awg"], "stranding": row["stranding"], "rating_A": row["current_rating_A"]}
        for conv in ("churchill_chu", "ieee738"):
            m = CableThermalModel(w.to_cable(20.0), convection_model=conv)
            out[f"rise_at_rating_{conv}"] = steady_T(m, row["current_rating_A"], 20.0) - 20.0
            out[f"I_40K_{conv}"] = m.ampacity(60.0, 20.0)
        rs.append(out)
    res["rs"] = rs
    # Alpha Wire PV wire (XLPE + PVC jacket, 90 °C), ratings at 30 °C
    pv = []
    for p in ("PV0830", "PV0630", "PV0430"):
        w = find_part(p)
        out = {"part": p, "awg": w.awg, "rating_A": w.rating_A}
        for conv in ("churchill_chu", "ieee738"):
            m = CableThermalModel(w.to_cable(30.0), convection_model=conv)
            out[f"T_at_rating_{conv}"] = steady_T(m, w.rating_A, 30.0)
            out[f"I_90C_{conv}"] = m.ampacity(90.0, 30.0)
        pv.append(out)
    res["alpha_pv"] = pv
    # General Cable Super Vu-Tron welding cable (90 °C, 30 °C ambient) - conductor diameter estimated
    meta = refs["general_cable_super_vutron"]
    svt = []
    for row in meta["rows"]:
        out = {"awg": row["awg"], "catalog": row["catalog"], "rating_A": row["current_rating_A"]}
        for pack in (0.75, 0.80, 0.85):
            w = from_supervutron(row, meta, pack)
            for conv in ("churchill_chu", "ieee738"):
                m = CableThermalModel(w.to_cable(30.0), convection_model=conv)
                if pack == 0.80:
                    out[f"T_at_rating_{conv}"] = steady_T(m, row["current_rating_A"], 30.0)
                out[f"I_90C_{conv}_pack{pack:.2f}"] = m.ampacity(90.0, 30.0)
        svt.append(out)
    res["super_vutron"] = svt
    return res


# ----------------------------------------------------------------------------- figures
def _awg_axis(ax, awgs):
    order = sorted(set(awgs), key=awg_sort_key, reverse=True)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    return {a: i for i, a in enumerate(order)}


def fig_coverage(rows, path):
    sim = [r for r in rows]
    groups = list(GROUPS) + ["not simulated"]
    awgs = sorted({r["awg"] for r in sim}, key=awg_sort_key, reverse=True)
    M = np.zeros((len(groups), len(awgs)), int)
    for r in sim:
        g = r["group"] if r["group"] in groups else "not simulated"
        M[groups.index(g), awgs.index(r["awg"])] += 1
    fig, ax = plt.subplots(figsize=(12, 4.6))
    ax.imshow(np.where(M > 0, np.log1p(M), np.nan), cmap="Blues", aspect="auto", vmin=0, vmax=np.log1p(M.max()) * 1.3)
    for i in range(len(groups)):
        for j in range(len(awgs)):
            if M[i, j]:
                ax.text(j, i, str(M[i, j]), ha="center", va="center", fontsize=8, color=INK)
    ax.set_xticks(range(len(awgs)))
    ax.set_xticklabels(awgs)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([f"{g} ({M[i].sum()})" for i, g in enumerate(groups)])
    ax.set_xlabel("Conductor size [AWG]")
    ax.grid(False)
    ax.set_title(f"Wire database: {len(rows)} Alpha Wire hook-up wire parts by insulation and AWG (count per cell)",
                 loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_ampacity(rows, path):
    sim = [r for r in rows if r.get("I_rated_A")]
    awgs = sorted({r["awg"] for r in sim}, key=awg_sort_key, reverse=True)
    pos = {a: i for i, a in enumerate(awgs)}
    fig, axs = plt.subplots(3, 3, figsize=(15, 11), sharex=True, sharey=True)
    for ax, (g, (_, col)) in zip(axs.flat, GROUPS.items()):
        for r in sim:
            ax.plot(pos[r["awg"]], r["I_105C_A"], "o", ms=3, color="#d8d7d2", zorder=1)
        mine = [r for r in sim if r["group"] == g]
        for r in mine:
            ax.plot(pos[r["awg"]], r["I_rated_A"], "o", ms=5.5, color=col, zorder=3)
        tms = sorted({int(r["T_max_C"]) for r in mine})
        ax.set_title(f"{g}: {len(mine)} parts, rated {', '.join(map(str, tms))} °C", loc="left", fontsize=9.5)
        ax.set_yscale("log")
    for ax in axs[-1]:
        ax.set_xticks(range(len(awgs)))
        ax.set_xticklabels(awgs, fontsize=7.5, rotation=90)
        ax.set_xlabel("AWG")
    for ax in axs[:, 0]:
        ax.set_ylabel("Current [A]")
    fig.suptitle(f"Model ampacity in still air at {AMBIENT:g} °C: coloured = conductor at the wire's rated temperature; "
                 "grey = all wires at 105 °C (Churchill-Chu convection, ε = 0.85)", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_material_effect(rows, path):
    """Same conductor temperature (105 °C): pure effect of insulation material and wall."""
    sim = [r for r in rows if r.get("I_105C_A")]
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    ax = axs[0]
    awgs = sorted({r["awg"] for r in sim}, key=awg_sort_key, reverse=True)
    pos = {a: i for i, a in enumerate(awgs)}
    for g, (_, col) in GROUPS.items():
        rr = [r for r in sim if r["group"] == g]
        if rr:
            ax.plot([pos[r["awg"]] for r in rr], [r["I_105C_A"] / r["metal_area_mm2"] for r in rr],
                    "o", ms=4.5, color=col, label=g, alpha=0.9)
    ax.set_xticks(range(len(awgs)))
    ax.set_xticklabels(awgs, fontsize=8, rotation=90)
    ax.set_xlabel("AWG")
    ax.set_ylabel("Current density at 105 °C conductor [A/mm²]")
    ax.set_yscale("log")
    ax.set_title("Allowable current density falls with size (surface/volume)", loc="left", fontsize=10)
    ax.legend(fontsize=7.5, ncol=2)
    ax = axs[1]
    for g, (_, col) in GROUPS.items():
        rr = [r for r in sim if r["group"] == g and r.get("tau_s") and np.isfinite(r["tau_s"])]
        if rr:
            ax.plot([r["OD_mm"] for r in rr], [r["tau_s"] for r in rr], "o", ms=4.5, color=col, label=g, alpha=0.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Outer diameter [mm]")
    ax.set_ylabel("Thermal time constant τ (63 %) [s]")
    ax.set_title("Heating time constant (step to the 40 K-rise current)", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_reference(res, path):
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.6))
    ax = axs[0]
    sets = [("Habia PTFE E (200 °C @ 40 °C amb.)", res["habia"]["rows"], "I_200C", VIOLET, "o"),
            ("RS PRO PTFE A (40 K rise)", res["rs"], "I_40K", AQUA, "s"),
            ("Alpha Wire PV (90 °C @ 30 °C)", res["alpha_pv"], "I_90C", YELLOW, "D"),
            ("Super Vu-Tron welding (90 °C @ 30 °C)", res["super_vutron"], "I_90C", GREY, "^")]
    lo, hi = 1.5, 600
    for lab, rows, key, col, mk in sets:
        x = [r["rating_A"] for r in rows]
        k_cc = key + "_churchill_chu" + ("_pack0.80" if "Vu-Tron" in lab else "")
        k_738 = key + "_ieee738" + ("_pack0.80" if "Vu-Tron" in lab else "")
        ax.plot(x, [r[k_cc] for r in rows], mk, color=col, ms=6.5, label=lab + " - Churchill-Chu")
        ax.plot(x, [r[k_738] for r in rows], mk, mfc="none", mec=col, ms=6.5, label=lab + " - Eq. 25")
    ax.plot([lo, hi], [lo, hi], color=INK2, lw=1)
    ax.fill_between([lo, hi], [lo * 0.8, hi * 0.8], [lo * 1.2, hi * 1.2], color="#e4e3df", alpha=0.6, lw=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Manufacturer current rating [A]")
    ax.set_ylabel("Model current for the same criterion [A]")
    ax.set_title("Model vs manufacturer ratings (band = ±20 %)", loc="left", fontsize=10)
    ax.legend(fontsize=7, loc="upper left")
    ax = axs[1]
    hb = res["habia"]["rows"]
    xs = range(len(hb))
    ax.plot(xs, [r["T_at_rating_churchill_chu"] for r in hb], "o-", color=VIOLET, label="Churchill-Chu")
    ax.plot(xs, [r["T_at_rating_ieee738"] for r in hb], "o--", mfc="none", color=VIOLET, label="Eq. 25 (paper)")
    ax.axhline(200, color=RED, lw=1)
    ax.text(0, 203, "PTFE rating 200 °C", color=RED, fontsize=8)
    ax.set_xticks(list(xs))
    ax.set_xticklabels([r["type"].replace("E ", "") for r in hb], rotation=90, fontsize=7.5)
    ax.set_xlabel("Habia type E (AWG + strands, e.g. 2219 = 22 AWG 19 strands)")
    ax.set_ylabel("Model conductor temperature at Habia rating [°C]")
    ax.set_title("Habia PTFE: steady temperature at the published current, 40 °C air", loc="left", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_derating(res, path):
    pub = res["habia"]["derating_published"]
    fig, ax = plt.subplots(figsize=(7, 4.4))
    T = sorted(int(k) for k in pub)
    ax.plot(T, [pub[str(t)] if str(t) in pub else pub[t] for t in T], "s-", color=INK, label="Habia published factor")
    for (t, d), col in zip(res["habia"]["derating_model"].items(), [VIOLET, AQUA, ORANGE]):
        ax.plot(T, [d[str(x)] if str(x) in d else d[x] for x in T], "o--", color=col, label=f"Model, {t} (200 °C conductor)")
    Ta = np.array(T, float)
    ax.plot(Ta, np.sqrt((200 - Ta) / 160), ":", color=GREY, label="√((200−T)/160)")
    ax.set_xlabel("Ambient air temperature [°C]")
    ax.set_ylabel("Current correction factor (40 °C = 1)")
    ax.set_title("Ambient correction: model vs Habia table", loc="left", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_transient_22awg(path):
    # all 22 AWG (7/30 or 19/34 strands), one representative per insulation family
    parts = [("3051", "PVC 0.016\" wall (UL1007)"), ("3071", "PVC 0.032\" wall (UL1015)"), ("5855/7", "PTFE 0.010\" wall"),
             ("39X2205", "Silicone rubber 0.020\" wall"), ("6713", "mPPE 0.0095\" wall (EcoWire)"),
             ("7024", "XLPE 0.016\" wall"), ("AZ221934", "ETFE 0.011\" wall"), ("2622", "Fluoroplastic 0.007\" (ThermoThin)")]
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    cols = [ORANGE, RED, VIOLET, BLUE, GREEN, YELLOW, AQUA, MAGENTA]
    out = {}
    for (p, lab), col in zip(parts, cols):
        try:
            w = find_part(p)
        except KeyError:
            continue
        m = CableThermalModel(w.to_cable(AMBIENT), convection_model=CONVECTION)
        r = m.simulate(8.0, ambient=AMBIENT, t_end=600.0, dt=0.5)
        axs[0].plot(r.t, r.T_center, color=col, label=f"{p} {w.awg} AWG {lab}, OD {w.outer_diameter * 1e3:.2f} mm")
        Tss = m.steady_state(8.0, AMBIENT)
        x = m.x / m.x[-1]
        axs[1].plot(m.x * 1e3, Tss, color=col)
        out[p] = {"awg": w.awg, "label": lab, "T_ss_C": float(Tss[0]), "T_surface_C": float(Tss[-1]),
                  "T_60s_C": float(np.interp(60, r.t, r.T_center))}
    axs[0].set_xlabel("Time [s]")
    axs[0].set_ylabel("Conductor temperature [°C]")
    axs[0].set_title(f"22 AWG wires, 8 A step, still air {AMBIENT:g} °C", loc="left", fontsize=10)
    axs[0].legend(fontsize=7.5)
    axs[1].set_xlabel("Radius [mm]")
    axs[1].set_ylabel("Steady temperature [°C]")
    axs[1].set_title("Steady radial profiles at 8 A (conductor | insulation)", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return out


def sensitivity(path):
    cases = [("2840/7", "32 AWG PTFE"), ("3051", "22 AWG PVC"), ("391450", "14 AWG silicone"),
             ("786133", "6 AWG PVC"), ("780001", "1/0 AWG PVC")]
    variants = [("k ×0.7", dict(k=0.7)), ("k ×1.3", dict(k=1.3)), ("ε = 0.70", dict(eps=0.70)), ("ε = 0.95", dict(eps=0.95)),
                ("Eq. 25 instead of Churchill-Chu", dict(conv="ieee738")), ("air 0.2 m/s", dict(v=0.2)),
                ("ambient 40 °C", dict(Ta=40.0))]
    res = {}
    fig, ax = plt.subplots(figsize=(10, 4.8))
    width = 0.11
    for i, (p, lab) in enumerate(cases):
        try:
            w = find_part(p)
        except KeyError:
            continue
        base_c = w.to_cable(AMBIENT)
        Tmax = float(w.temp_max_C)
        I0 = CableThermalModel(base_c, convection_model=CONVECTION).ampacity(Tmax, AMBIENT)
        res[p] = {"label": lab, "I0": I0}
        for j, (vl, kw) in enumerate(variants):
            c = base_c
            if "k" in kw:
                from .properties import Material
                scale = kw["k"]
                c = c.with_(insulation_material=Material(c.insulation_material.name, c.insulation_material.density,
                                                         c.insulation_material.cp, lambda T, f=c.insulation_material.k: scale * f(T)),
                            jacket_material=Material(c.jacket_material.name, c.jacket_material.density,
                                                     c.jacket_material.cp, lambda T, f=c.jacket_material.k: scale * f(T)))
            if "eps" in kw:
                c = c.with_(emissivity=kw["eps"])
            m = CableThermalModel(c, convection_model=kw.get("conv", CONVECTION), air_speed=kw.get("v", 0.0))
            I = m.ampacity(Tmax, kw.get("Ta", AMBIENT))
            res[p][vl] = I / I0 - 1
            ax.bar(i + (j - 3) * width, 100 * (I / I0 - 1), width, color=[ORANGE, YELLOW, BLUE, VIOLET, RED, AQUA, GREY][j],
                   label=vl if i == 0 else None)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels([f"{p}\n{lab}" for p, lab in cases], fontsize=8.5)
    ax.set_ylabel("Change of ampacity [%]")
    ax.set_title("Ampacity sensitivity (rated conductor temperature, 30 °C air)", loc="left", fontsize=10)
    ax.legend(fontsize=7.5, ncol=4)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return res


# ----------------------------------------------------------------------------- report
def write_report(rows, cmp, tr, sens, path):
    L = []
    w = L.append
    sim = [r for r in rows if r.get("I_rated_A")]
    w("# Wire database and test campaign\n")
    w("Generated by `python -m cable_thermal wire-study`. Model: the paper's radial FDM/TDMA model "
      "(9 nodes per material, still air, emissivity 0.85, copper α = 0.00393 1/K) with Churchill–Chu natural "
      f"convection (Eq. 25 results in brackets where shown). Ambient {AMBIENT:g} °C unless stated. "
      "Data sources and modelling choices: docs/WIRE_DATABASE.md.\n")
    w("## 1. Database coverage\n")
    w(f"{len(rows)} Alpha Wire parts ({min(rows, key=lambda r: -awg_sort_key(r['awg']))['awg']} AWG to "
      f"{min(rows, key=lambda r: awg_sort_key(r['awg']))['awg']} AWG), {len(sim)} simulated "
      f"({len(rows) - len(sim)} mica-glass/fiberglass parts have no thermal data). Plus 25 Habia, 3 RS PRO and 8 General Cable "
      "rows with published current ratings.\n")
    c = Counter(r["group"] for r in rows)
    w("| Insulation group | parts | AWG range | rated temperatures [°C] |")
    w("|---|---|---|---|")
    for g in list(GROUPS) + ["not simulated"]:
        rr = [r for r in rows if r["group"] == g]
        if rr:
            a = sorted({r["awg"] for r in rr}, key=awg_sort_key)
            tm = sorted({int(r["T_max_C"]) for r in rr if r["T_max_C"]})
            w(f"| {g} | {len(rr)} | {a[-1]} – {a[0]} | {', '.join(map(str, tm))} |")
    w("")
    w("## 2. Model ampacity by AWG (median over parts, rated conductor temperature, 30 °C still air)\n")
    awgs = sorted({r["awg"] for r in sim}, key=awg_sort_key, reverse=True)
    groups = [g for g in GROUPS if any(r["group"] == g for r in sim)]
    w("| AWG | " + " | ".join(groups) + " |")
    w("|---" * (len(groups) + 1) + "|")
    for a in awgs:
        cells = []
        for g in groups:
            v = [r["I_rated_A"] for r in sim if r["awg"] == a and r["group"] == g]
            cells.append(f"{np.median(v):.3g}" if v else "")
        w(f"| {a} | " + " | ".join(cells) + " |")
    w("\nAll per-part results (ampacity at rated temperature with both convection laws, at 105 °C, for a 40 K rise, "
      "time constant, insulation temperature drop) are in `wire_results.csv`.\n")
    w("## 3. Comparison with manufacturer current ratings\n")
    w("### Habia PTFE type E (NEMA HP-3, 200 °C), 'Amps @ 40 °C'\n")
    w("| Type | AWG | rating [A] | model T at rating [°C] CC (Eq.25) | model I for 200 °C [A] CC (Eq.25) | rating / model |")
    w("|---|---|---|---|---|---|")
    for r in cmp["habia"]["rows"]:
        w(f"| {r['type']} | {r['awg']} | {r['rating_A']} | {r['T_at_rating_churchill_chu']:.0f} ({r['T_at_rating_ieee738']:.0f}) | "
          f"{r['I_200C_churchill_chu']:.1f} ({r['I_200C_ieee738']:.1f}) | {r['rating_A'] / r['I_200C_churchill_chu']:.2f} |")
    ratios = [r["rating_A"] / r["I_200C_churchill_chu"] for r in cmp["habia"]["rows"]]
    Tr = [r["T_at_rating_churchill_chu"] for r in cmp["habia"]["rows"]]
    w(f"\nAt the published currents the model gives conductor temperatures of {min(Tr):.0f}–{max(Tr):.0f} °C "
      f"(median {np.median(Tr):.0f} °C) in 40 °C air: the Habia ratings correspond to a conductor close to the 200 °C "
      f"insulation limit (rating / model-200 °C current = {min(ratios):.2f}–{max(ratios):.2f}, median {np.median(ratios):.2f}).\n")
    w("Ambient correction factors (current at 200 °C conductor relative to 40 °C air):\n")
    pub = cmp["habia"]["derating_published"]
    T = sorted(int(k) for k in pub)
    w("| T air [°C] | " + " | ".join(map(str, T)) + " |")
    w("|---" * (len(T) + 1) + "|")
    w("| Habia | " + " | ".join(f"{pub[str(t)] if str(t) in pub else pub[t]:.2f}" for t in T) + " |")
    for t, d in cmp["habia"]["derating_model"].items():
        w(f"| model {t} | " + " | ".join(f"{d[str(x)] if str(x) in d else d[x]:.2f}" for x in T) + " |")
    w("")
    w("### RS PRO PTFE type A (BS 3G 210), rating = 40 K rise, single wire in free air (20 °C air assumed)\n")
    w("| AWG | stranding | rating [A] | model rise at rating [K] CC (Eq.25) | model current for 40 K [A] CC (Eq.25) |")
    w("|---|---|---|---|---|")
    for r in cmp["rs"]:
        w(f"| {r['awg']} | {r['stranding']} | {r['rating_A']:g} | {r['rise_at_rating_churchill_chu']:.1f} ({r['rise_at_rating_ieee738']:.1f}) | "
          f"{r['I_40K_churchill_chu']:.1f} ({r['I_40K_ieee738']:.1f}) |")
    w("")
    w("### Alpha Wire PV wire (XLPE + PVC jacket, 90 °C), datasheet current at 30 °C\n")
    w("| Part | AWG | rating [A] | model T at rating [°C] CC (Eq.25) | model I for 90 °C [A] CC (Eq.25) |")
    w("|---|---|---|---|---|")
    for r in cmp["alpha_pv"]:
        w(f"| {r['part']} | {r['awg']} | {r['rating_A']} | {r['T_at_rating_churchill_chu']:.1f} ({r['T_at_rating_ieee738']:.1f}) | "
          f"{r['I_90C_churchill_chu']:.0f} ({r['I_90C_ieee738']:.0f}) |")
    w("")
    w("### General Cable Super Vu-Tron welding cable (90 °C, NEC 400.5(B) basis, 30 °C) – 6 AWG to 4/0 AWG\n")
    w("Conductor diameter not published: estimated from the metal area with a packing factor 0.80 (range 0.75–0.85 shown); "
      "rubber modelled with EPR data.\n")
    w("| AWG | rating [A] | model T at rating [°C] CC (Eq.25) | model I for 90 °C [A], packing 0.75 / 0.80 / 0.85 (CC) |")
    w("|---|---|---|---|")
    for r in cmp["super_vutron"]:
        w(f"| {r['awg']} | {r['rating_A']} | {r['T_at_rating_churchill_chu']:.1f} ({r['T_at_rating_ieee738']:.1f}) | "
          + " / ".join(f"{r[f'I_90C_churchill_chu_pack{p:.2f}']:.0f}" for p in (0.75, 0.80, 0.85)) + " |")
    w("")
    w("## 4. Same conductor (22 AWG), different insulation: 8 A step in 30 °C still air\n")
    w("| Part | insulation | steady conductor [°C] | surface [°C] | conductor after 60 s [°C] |")
    w("|---|---|---|---|---|")
    for p, v in tr.items():
        w(f"| {p} | {v['label']} | {v['T_ss_C']:.1f} | {v['T_surface_C']:.1f} | {v['T_60s_C']:.1f} |")
    w("")
    w("## 5. Sensitivity of the ampacity\n")
    vl = [k for k in next(iter(sens.values())) if k not in ("label", "I0")]
    w("| Wire | ampacity [A] | " + " | ".join(vl) + " |")
    w("|---" * (len(vl) + 2) + "|")
    for p, v in sens.items():
        w(f"| {p} {v['label']} | {v['I0']:.1f} | " + " | ".join(f"{100 * v[k]:+.1f} %" for k in vl) + " |")
    w("")
    w("## 6. Data-quality flags\n")
    flagged = [r for r in rows if r["notes"]]
    for r in flagged:
        w(f"* {r['part']} ({r['awg']} AWG {r['insulation']}): {r['notes']}")
    rho = [(r["part"], r["awg"], r["R20_ohm_per_km"] * 1e-3 * r["metal_area_mm2"] * 1e-6) for r in rows
           if r["R20_source"] == "datasheet DCR"]
    low = [x for x in rho if x[2] < 0.97 * 1.7241e-8]   # more than 3 % below annealed copper
    for p, a, rr in low:
        w(f"* {p} ({a} AWG): datasheet DCR × metal area of the stated stranding = {rr * 1e8:.2f}e-8 Ω·m, "
          "below the resistivity of annealed copper (1.724e-8) – the DCR or the stranding on the datasheet is inconsistent; "
          "the datasheet DCR is used.")
    w(f"* Median implied resistivity of all parts with a DCR: {np.median([x[2] for x in rho]) * 1e8:.2f}e-8 Ω·m "
      "(tinned/stranded conductors, nominal values).")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


def run_wire_study(out="results/wires"):
    _style()
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    F = lambda n: os.path.join(out, "figures", n)
    rows = run_campaign(os.path.join(out, "wire_results.csv"))
    cmp = run_comparisons()
    fig_coverage(rows, F("fig20_wire_database_coverage.png"))
    fig_ampacity(rows, F("fig21_wire_ampacity_by_material.png"))
    fig_material_effect(rows, F("fig22_current_density_time_constant.png"))
    fig_reference(cmp, F("fig23_model_vs_manufacturer_ratings.png"))
    fig_derating(cmp, F("fig24_habia_ambient_correction.png"))
    tr = fig_transient_22awg(F("fig25_22awg_insulation_transient.png"))
    sens = sensitivity(F("fig26_wire_ampacity_sensitivity.png"))
    with open(os.path.join(out, "comparisons.json"), "w") as fh:
        json.dump({"comparisons": cmp, "transient_22awg": tr, "sensitivity": sens}, fh, indent=1, default=float)
    write_report(rows, cmp, tr, sens, os.path.join(out, "RESULTS_WIRES.md"))
    return rows, cmp
