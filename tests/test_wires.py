"""Tests of the wire database, its conversion to model inputs and the thin-wire convection law."""
import math

import numpy as np
import pytest

from cable_thermal import CableThermalModel
from cable_thermal.heat_transfer import h_churchill_chu, h_natural, nusselt_churchill_chu, rayleigh_prandtl
from cable_thermal.properties import constant_material
from cable_thermal.reference import analytic_steady_state_constant_props
from cable_thermal.wires import (IN, KFT, all_alpha_wires, alpha_wire_parts, awg_diameter_m, find_part, from_rs,
                                 material_db, parse_stranding, reference_ratings)


def test_awg_diameters_astm_b258():
    assert np.isclose(awg_diameter_m(36), 0.127e-3)
    assert np.isclose(awg_diameter_m(10), 2.588e-3, rtol=1e-3)
    assert np.isclose(awg_diameter_m(30), 0.2546e-3, rtol=1e-3)
    assert np.isclose(awg_diameter_m("1/0"), 8.252e-3, rtol=1e-3)
    assert np.isclose(awg_diameter_m("4/0"), 11.684e-3, rtol=1e-3)


@pytest.mark.parametrize("s,awg,n,d", [("7/30", "22", 7, awg_diameter_m(30)), ("SOLID", "22", 1, awg_diameter_m(22)),
                                       ("133(7x19)/29", "8", 133, awg_diameter_m(29)),
                                       ("259(37x7)25", "1", 259, awg_diameter_m(25)),
                                       ("19/.0591", "2", 19, 0.0591 * IN), ("665 (7x95)/30", "2", 665, awg_diameter_m(30))])
def test_parse_stranding(s, awg, n, d):
    nn, dd = parse_stranding(s, awg)
    assert nn == n and np.isclose(dd, d)


def test_database_contents():
    parts = alpha_wire_parts()
    assert len(parts) == 337
    assert not any("parse_error" in p for p in parts)
    awgs = {p["awg"] for p in parts}
    assert {"36", "30", "22", "10", "4", "2", "1", "1/0"} <= awgs
    mats = {l["material"] for p in parts for l in p["insulation_layers"]}
    assert {"PVC", "Silicone Rubber", "PTFE", "XLPE", "MPPE", "ETFE", "PVDF"} <= mats
    mapping = material_db()["alpha_wire_mapping"]
    assert mats <= set(mapping)                                    # every material is mapped
    for key in {v for v in mapping.values() if v}:
        m = material_db()["materials"][key]
        for src in ("k_source", "density_source", "cp_source"):
            assert m[src] in material_db()["_sources"] or m[src] == "paper" or m[src].startswith("assumed")


def test_alpha_record_conversion():
    w = find_part("3051")                                          # 22 AWG 7/30 tinned Cu, 0.016" PVC, OD 0.062"
    assert w.awg == "22" and w.layers[0][0] == "PVC"
    assert np.isclose(w.conductor_diameter, 0.030 * IN) and np.isclose(w.outer_diameter, 0.062 * IN)
    assert np.isclose(w.r20, 16.2 / KFT)                            # datasheet DCR 16.2 ohm/1000 ft
    assert np.isclose(w.metal_area, 7 * math.pi / 4 * awg_diameter_m(30) ** 2)
    c = w.to_cable(30.0)
    assert np.isclose(c.outer_diameter, w.outer_diameter)
    assert np.isclose(c.resistance_per_length(20.0), w.r20)
    assert c.insulation_material is c.jacket_material or c.insulation_material.name == c.jacket_material.name
    # every simulatable record yields a geometry whose OD is the datasheet OD
    for w in all_alpha_wires():
        if w.simulatable:
            c = w.to_cable()
            assert np.isclose(c.outer_diameter, w.outer_diameter, rtol=1e-9)
            assert c.insulation_thickness > 0 and c.jacket_thickness > 0


def test_single_insulation_split_matches_analytic():
    """A single insulation layer represented as two equal sub-layers = exact cylinder conduction."""
    w = find_part("3071")
    c = w.to_cable(30.0).with_(insulation_material=constant_material("i", 1400, 1000, 0.2),
                               jacket_material=constant_material("i", 1400, 1000, 0.2),
                               conductor_material=constant_material("cu", 8960, 390, 398), alpha=0.0)
    m = CableThermalModel(c, nodes_per_material=(40, 40, 40))
    T = m.steady_state(10.0)
    q = 100.0 * c.resistance_per_length(20.0)
    T0, TAB, _ = analytic_steady_state_constant_props(c, q, T[-1], 398, 0.2, 0.2)
    assert abs(T[m.probes["conductor_insulation"]] - TAB) < 0.002 * (T0 - T[-1]) + 1e-6


def test_churchill_chu_formula_and_limits():
    # hand evaluation of Churchill & Chu at Ra = 1e4, Pr = 0.7
    nu = (0.6 + 0.387 * 1e4 ** (1 / 6) / (1 + (0.559 / 0.7) ** (9 / 16)) ** (8 / 27)) ** 2
    assert np.isclose(nusselt_churchill_chu(1e4, 0.7), nu) and np.isclose(nu, 4.366, atol=2e-3)
    assert np.isclose(nusselt_churchill_chu(0.0, 0.7), 0.36)        # conduction limit Ra -> 0
    Ra, Pr, k = rayleigh_prandtl(60.0, 20.0, 0.0147)
    assert 0.69 < Pr < 0.72 and 0.026 < k < 0.029
    # thin wires: power-law Eq. 25 under-predicts h; for the paper's cables both agree within ~15 %
    assert h_churchill_chu(60, 20, 0.5e-3) > 1.6 * h_natural(60, 20, 0.5e-3, model="ieee738")
    assert abs(h_churchill_chu(60, 20, 0.0147) / h_natural(60, 20, 0.0147, model="ieee738") - 1) < 0.15


def test_rs_ptfe_ratings_correspond_to_about_40K():
    meta = reference_ratings()["rs_ptfe_type_a"]
    for row in meta["rows"]:
        m = CableThermalModel(from_rs(row, meta).to_cable(20.0), convection_model="churchill_chu")
        rise = m.steady_state(row["current_rating_A"], 20.0)[0] - 20.0
        assert 30.0 < rise < 50.0


def test_wire_ampacity_regression():
    m = CableThermalModel(find_part("3051").to_cable(30.0), convection_model="churchill_chu")
    assert 11.0 < m.ampacity(105.0, 30.0) < 12.2
    m = CableThermalModel(find_part("PV0430").to_cable(30.0), convection_model="churchill_chu")
    assert 140.0 < m.ampacity(90.0, 30.0) < 175.0                  # datasheet: 140 A at 30 °C
