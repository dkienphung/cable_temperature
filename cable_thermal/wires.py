"""Wire database: manufacturer hook-up wire / cable data turned into model-ready CableSpec objects.

Data (``cable_thermal/data``):
* ``wires/alphawire_hookup.json``   - 337 Alpha Wire hook-up wire parts, 36 AWG ... 1/0 AWG
* ``wires/reference_ratings.json``  - third-party data with published current ratings
                                      (Habia PTFE type E, RS PRO PTFE type A, General Cable Super Vu-Tron)
* ``insulation_materials.json``     - k, density, c_p of the insulation materials, with sources

How a datasheet becomes a model (all choices documented in docs/WIRE_DATABASE.md):
* conductor: datasheet diameter D1; metal area from the stranding (ASTM B258 AWG diameters);
  r20 = datasheet DC resistance (ohm/1000 ft at 20 °C) when given, alpha = 0.00393 1/K
  (annealed copper, 20 °C); mass = 8960 kg/m3 x metal area (plating neglected);
* insulation: one layer of wall t is represented as two equal sub-layers of the same material
  (the paper's 3-material scheme with B = C is exactly a single layer); insulation + jacket map
  to B and C. The outermost wall is adjusted so that the model's outer diameter equals the
  datasheet OD (the surface that convects and radiates);
* surface: still air, emissivity 0.85 (paper), Churchill & Chu natural convection by default
  (Eq. 25 of the paper is a power law outside its range for thin wires).
"""
from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from typing import Dict, List, Optional

import numpy as np

from .cables import CableSpec
from .properties import COPPER_DENSITY, Material, constant_material, copper, pvc, xlpe

IN = 25.4e-3                      # m per inch
KFT = 304.8                       # m per 1000 ft
RHO_CU_20 = 1.7241e-8             # ohm m, annealed copper (IEC 60028 / ASTM B3)
ALPHA_CU_20 = 0.00393             # 1/K, annealed copper at 20 °C
DEFAULT_EMISSIVITY = 0.85         # paper


def awg_diameter_m(awg) -> float:
    """ASTM B258 definition: d = 0.127 mm x 92^((36 - n)/39); n = -1, -2, ... for 1/0, 2/0, ..."""
    s = str(awg).strip().upper()
    if s.endswith("/0"):
        n = 1 - int(s.split("/")[0])
    else:
        n = int(s)
    return 0.127e-3 * 92 ** ((36 - n) / 39)


def awg_area_m2(awg) -> float:
    return math.pi / 4 * awg_diameter_m(awg) ** 2


def parse_stranding(stranding: str, awg: str):
    """Return (number of strands, strand diameter [m]) from strings such as '7/30', 'SOLID',
    '133(7x19)/29', '259(37x7)25', '19/.0591', '665 (7x95)/30'."""
    s = stranding.strip().upper().replace(" ", "")
    if s == "SOLID":
        return 1, awg_diameter_m(awg)
    m = re.match(r"^(\d+)(?:\([^)]*\))?/?(\.?\d+(?:\.\d+)?)$", s)
    if not m:
        raise ValueError(f"cannot parse stranding {stranding!r}")
    n, g = int(m.group(1)), m.group(2)
    d = float(g) * IN if g.startswith(".") else awg_diameter_m(int(g))
    return n, d


@lru_cache(maxsize=None)
def _load(name: str):
    if name.startswith("wires/"):
        return json.loads(resources.files("cable_thermal.data.wires").joinpath(name[6:]).read_text())
    return json.loads(resources.files("cable_thermal.data").joinpath(name).read_text())


def material_db() -> dict:
    return _load("insulation_materials.json")


def insulation_material(key: str) -> Material:
    """Material object for a key of insulation_materials.json."""
    m = material_db()["materials"][key]
    if m["model"] == "paper_curves":
        return {"PVC": pvc, "XLPE": xlpe}[key](m["density"])
    return constant_material(key, m["density"], m["cp"], m["k"])


def map_alpha_material(name: str) -> Optional[str]:
    mp = material_db()["alpha_wire_mapping"]
    if name not in mp:
        raise KeyError(f"no material mapping for {name!r}")
    return mp[name]


def alpha_wire_parts() -> List[dict]:
    return _load("wires/alphawire_hookup.json")["parts"]


def reference_ratings() -> dict:
    return _load("wires/reference_ratings.json")


@dataclass
class WireRecord:
    """Normalised wire description (SI units) shared by all data sources."""
    source: str
    part_number: str
    awg: str
    conductor_diameter: float          # m
    metal_area: float                  # m²
    r20: float                         # ohm/m at 20 °C
    r20_source: str
    layers: list                       # [(material key, thickness m)], inside -> outside
    outer_diameter: float              # m
    temp_max_C: Optional[float]
    insulation_label: str
    notes: tuple = ()
    rating_A: Optional[float] = None
    rating_basis: Optional[str] = None
    raw: Optional[dict] = None

    @property
    def simulatable(self) -> bool:
        return all(k is not None for k, _ in self.layers)

    def to_cable(self, ambient: float = 30.0, emissivity: float = DEFAULT_EMISSIVITY,
                 elevation: float = 0.0) -> CableSpec:
        if not self.simulatable:
            raise ValueError(f"{self.part_number}: no thermal data for {self.insulation_label}")
        if len(self.layers) == 1:
            key, t = self.layers[0]
            mats = [insulation_material(key)] * 2
            th = [t / 2, t / 2]
        elif len(self.layers) == 2:
            mats = [insulation_material(k) for k, _ in self.layers]
            th = [t for _, t in self.layers]
        else:
            raise ValueError("at most two non-conductor layers are supported")
        return CableSpec(
            name=self.part_number, designation=f"{self.source} {self.part_number} {self.awg} AWG {self.insulation_label}",
            conductor_diameter=self.conductor_diameter, insulation_thickness=th[0], jacket_thickness=th[1],
            effective_cross_section=self.metal_area, resistivity_T0=self.r20 * self.metal_area, alpha=ALPHA_CU_20,
            T0=20.0, conductor_mass_per_length=COPPER_DENSITY * self.metal_area, conductor_material=copper(),
            insulation_material=mats[0], jacket_material=mats[1], emissivity=emissivity,
            ambient_temperature=ambient, elevation=elevation,
            max_continuous_temperature=float(self.temp_max_C) if self.temp_max_C is not None else 90.0,
            notes=self.notes)


def _fit_layers(d: float, od: float, layers: list):
    """Make the outer diameter equal to the datasheet OD by adjusting the outermost wall."""
    inner = sum(t for _, t in layers[:-1])
    t_last = (od - d) / 2 - inner
    nominal = layers[-1][1]
    if t_last <= 0:
        return layers, abs(d + 2 * sum(t for _, t in layers) - od) / od
    out = layers[:-1] + [(layers[-1][0], t_last)]
    return out, abs(t_last - nominal) / max(nominal, 1e-12)


def from_alpha(rec: dict) -> WireRecord:
    awg = rec["awg"]
    n, ds = parse_stranding(rec["stranding"], awg)
    area = n * math.pi / 4 * ds**2
    d = rec["conductor_diameter_in"] * IN
    notes = []
    if rec.get("dcr_ohm_per_kft_20C"):
        r20, r_src = rec["dcr_ohm_per_kft_20C"] / KFT, "datasheet DCR"
    else:
        r20, r_src = RHO_CU_20 / area, "computed (no DCR on datasheet)"
        notes.append("DC resistance computed from area with annealed-copper resistivity")
    layers = [(map_alpha_material(l["material"]), l["wall_in"] * IN) for l in rec.get("insulation_layers", [])]
    label = " / ".join(l["material"] for l in rec.get("insulation_layers", []))
    od = rec.get("insulation_od_in")
    if rec.get("jacket"):
        layers.append((map_alpha_material(rec["jacket"]["material"]), rec["jacket"]["wall_in"] * IN))
        label += " + " + rec["jacket"]["material"] + " jacket"
        od = rec["jacket"]["od_in"]
    od = od * IN if od else d + 2 * sum(t for _, t in layers)
    layers, dev = _fit_layers(d, od, layers)
    if dev > 0.05:
        notes.append(f"outer wall adjusted by {100 * dev:.0f}% to match the datasheet OD")
    return WireRecord("Alpha Wire", rec["part_number"], awg, d, area, r20, r_src, layers, od, rec.get("temp_max_C"),
                      label, tuple(notes), rec.get("current_rating_A"),
                      f"datasheet: {rec.get('current_rating_A')} A per conductor @ {rec.get('current_rating_ambient_C')} °C"
                      if rec.get("current_rating_A") else None, rec)


def from_habia(row: dict, meta: dict) -> WireRecord:
    d, od = row["conductor_diameter_mm"] * 1e-3, row["od_mm"] * 1e-3
    area = row["strands"] * math.pi / 4 * (row["strand_diameter_mm"] * 1e-3) ** 2
    return WireRecord("Habia", row["type"], row["awg"], d, area, row["dcr_ohm_per_km_20C"] * 1e-3, "datasheet DCR",
                      [("PTFE", (od - d) / 2)], od, meta["temp_max_C"], "PTFE", (), row["current_rating_A"],
                      meta["rating_basis"], row)


def from_rs(row: dict, meta: dict) -> WireRecord:
    d = row["conductor_diameter_mm"] * 1e-3
    od = 0.5 * (row["od_mm_min"] + row["od_mm_max"]) * 1e-3
    n, dsm = row["stranding"].split("/")
    area = int(n) * math.pi / 4 * (float(dsm) * 1e-3) ** 2
    return WireRecord("RS PRO", row["product"].split("(Type A)")[-1].strip(), row["awg"], d, area,
                      row["dcr_ohm_per_km_20C_max"] * 1e-3, "datasheet DCR (maximum)", [("PTFE", (od - d) / 2)], od,
                      meta["temp_max_C"], "PTFE", ("OD = mid-point of the published min/max",), row["current_rating_A"],
                      meta["rating_basis"], row)


def from_supervutron(row: dict, meta: dict, packing: float = 0.80) -> WireRecord:
    """Welding cable: the sheet gives OD and stranding only. Conductor diameter estimated from the metal
    area with a packing factor (fine-stranded rope ~0.75-0.85); resistance from annealed copper."""
    n, ds = row["strands"], awg_diameter_m(row["strand_awg"])
    area = n * math.pi / 4 * ds**2
    d = math.sqrt(4 * area / (math.pi * packing))
    od = row["od_in"] * IN
    return WireRecord("General Cable", "Super Vu-Tron " + row["catalog"], row["awg"], d, area, 1.02 * RHO_CU_20 / area,
                      "computed: annealed Cu x 1.02 stranding factor (not on datasheet)", [("EPR", (od - d) / 2)], od,
                      meta["temp_max_C"], "rubber (modelled as EPR)",
                      (f"conductor diameter estimated with packing factor {packing}", "insulation modelled with EPR (IEC) data"),
                      row["current_rating_A"], meta["rating_basis"], row)


def all_alpha_wires() -> List[WireRecord]:
    return [from_alpha(r) for r in alpha_wire_parts()]


def find_part(part: str) -> WireRecord:
    p = part.strip().upper()
    for r in alpha_wire_parts():
        if r["part_number"].upper() == p:
            return from_alpha(r)
    refs = reference_ratings()
    for row in refs["habia_E_ptfe"]["rows"]:
        if row["type"].replace(" ", "").upper() in (p.replace(" ", ""), "HABIA" + p.replace(" ", "")):
            return from_habia(row, refs["habia_E_ptfe"])
    raise KeyError(f"part {part!r} not found")


def awg_sort_key(awg: str) -> float:
    s = str(awg)
    return 1 - int(s.split("/")[0]) if s.endswith("/0") else int(s)
