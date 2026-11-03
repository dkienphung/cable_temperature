"""Temperature-dependent material properties used by the cable model.

The paper (Riba & Llauradó, Materials 2022, 15, 6814) gives the temperature
dependence of the specific heat capacity (Fig. 2) and thermal conductivity
(Fig. 3) of copper, XLPE and PVC only as plots between 15/20 °C and 75 °C.
Those curves were digitized (see ``data/material_properties_digitized.json``
and ``tools/``) and are interpolated linearly here.

Outside the digitized range the value is held constant at the nearest end point
("clamp"), which is the most conservative choice that does not invent data.
This matters only for scenarios that leave 20–75 °C (rated temperature 90 °C,
short circuits up to 250 °C) and is flagged in the documentation.

Data that the paper does NOT give (mass densities of XLPE and PVC, density of
copper for the 150 mm² cable) are explicit, documented assumptions in
:mod:`cable_thermal.cables`.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from importlib import resources
from typing import Callable, Sequence

import numpy as np

_DATA = json.loads(
    resources.files("cable_thermal.data").joinpath("material_properties_digitized.json").read_text()
)


@dataclass(frozen=True)
class Tabulated:
    """Piecewise-linear property table T [°C] -> value, clamped outside the range."""

    T: tuple
    values: tuple
    name: str = ""

    def __call__(self, T_celsius):
        return np.interp(T_celsius, self.T, self.values)

    @classmethod
    def constant(cls, value: float, name: str = "") -> "Tabulated":
        return cls((0.0,), (float(value),), name)

    @classmethod
    def from_pairs(cls, pairs: Sequence[Sequence[float]], name: str = "") -> "Tabulated":
        arr = np.asarray(pairs, dtype=float)
        order = np.argsort(arr[:, 0])
        return cls(tuple(arr[order, 0]), tuple(arr[order, 1]), name)

    @property
    def T_range(self):
        return (self.T[0], self.T[-1])


def digitized(key: str) -> Tabulated:
    """Property curve digitized from the paper (keys: copper_cp, xlpe_cp, pvc_cp,
    copper_k, xlpe_k, pvc_k)."""
    return Tabulated.from_pairs(_DATA[key]["data"], name=key)


@dataclass(frozen=True)
class Material:
    """Thermal material: density [kg/m³], specific heat c_p(T) [J/(kg K)],
    thermal conductivity k(T) [W/(m K)]."""

    name: str
    density: float
    cp: Callable = field(repr=False)
    k: Callable = field(repr=False)

    def rho_cp(self, T):
        return self.density * self.cp(T)


# --- materials of the paper -------------------------------------------------
COPPER_DENSITY = 8960.0      # kg/m³ (not given in the paper; standard value for Cu)
XLPE_DENSITY = 920.0         # kg/m³ (NOT given in the paper; typical XLPE 0.92 g/cm³)
PVC_DENSITY = 1400.0         # kg/m³ (NOT given in the paper; typical PVC 1.38-1.45 g/cm³)


def copper(density: float = COPPER_DENSITY) -> Material:
    return Material("copper", density, digitized("copper_cp"), digitized("copper_k"))


def xlpe(density: float = XLPE_DENSITY) -> Material:
    return Material("XLPE", density, digitized("xlpe_cp"), digitized("xlpe_k"))


def pvc(density: float = PVC_DENSITY) -> Material:
    return Material("PVC", density, digitized("pvc_cp"), digitized("pvc_k"))


def constant_material(name: str, density: float, cp: float, k: float) -> Material:
    """Constant-property material (used for verification against analytical solutions)."""
    return Material(name, density, Tabulated.constant(cp, name + "_cp"), Tabulated.constant(k, name + "_k"))


def material_at_reference(m: Material, T_ref: float = 20.0) -> Material:
    """Freeze a temperature-dependent material at T_ref (sensitivity studies)."""
    return constant_material(m.name + f"@{T_ref:g}C", m.density, float(m.cp(T_ref)), float(m.k(T_ref)))
