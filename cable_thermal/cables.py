"""Cable data of the paper (Tables 1-4) and helpers to define new cables.

Geometry (paper Fig. 1): D1 = conductor diameter, D2 = D1 + 2 t_ins,
D3 = D2 + 2 t_jacket. R_A = D1/2, R_B = D2/2, R_C = D3/2.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import List, Optional, Tuple

import numpy as np

from .properties import Material, copper, pvc, xlpe


@dataclass(frozen=True)
class CableSpec:
    name: str
    designation: str
    conductor_diameter: float          # D1 [m]  ("outer conductor diameter")
    insulation_thickness: float        # [m]
    jacket_thickness: float            # [m]
    effective_cross_section: float     # S_eff [m²]  (electrical cross-section)
    resistivity_T0: float              # rho_e(T0) [Ohm m]
    alpha: float                       # temperature coefficient of resistance [1/K]
    T0: float = 20.0                   # reference temperature of rho_e [°C]
    conductor_mass_per_length: Optional[float] = None   # m [kg/m]; None -> rho_Cu * S_eff
    conductor_material: Material = field(default_factory=copper)
    insulation_material: Material = field(default_factory=xlpe)
    jacket_material: Material = field(default_factory=pvc)
    emissivity: float = 0.85           # paper, Section 2.4
    ambient_temperature: float = 20.0  # [°C]
    test_steps: Tuple[Tuple[float, float], ...] = ()   # ((I_rms [A], duration [s]), ...)
    max_continuous_temperature: float = 90.0
    short_circuit_temperature: float = 250.0
    notes: Tuple[str, ...] = ()

    # --- geometry ---------------------------------------------------------
    @property
    def R_A(self):
        return 0.5 * self.conductor_diameter

    @property
    def R_B(self):
        return self.R_A + self.insulation_thickness

    @property
    def R_C(self):
        return self.R_B + self.jacket_thickness

    @property
    def outer_diameter(self):
        return 2.0 * self.R_C

    @property
    def conductor_geometric_area(self):
        """pi R_A^2 - area of the conductor region of the radial model."""
        return np.pi * self.R_A**2

    @property
    def conductor_mass(self) -> float:
        if self.conductor_mass_per_length is not None:
            return self.conductor_mass_per_length
        return self.conductor_material.density * self.effective_cross_section

    @property
    def conductor_density_effective(self) -> float:
        """Homogenized density of the conductor region so that
        sum(rho_A S_i,A) over the conductor equals the conductor mass m [kg/m]."""
        return self.conductor_mass / self.conductor_geometric_area

    # --- electrical ---------------------------------------------------------
    def resistance_per_length(self, T):
        """Eq. (2): r(T) = r_T0 [1 + alpha (T - T0)],  r_T0 = rho_e(T0) / S_eff  [Ohm/m]."""
        r0 = self.resistivity_T0 / self.effective_cross_section
        return r0 * (1.0 + self.alpha * (np.asarray(T, dtype=float) - self.T0))

    def with_(self, **kw) -> "CableSpec":
        return replace(self, **kw)

    def summary(self) -> dict:
        return {
            "name": self.name,
            "designation": self.designation,
            "D1_conductor_mm": self.conductor_diameter * 1e3,
            "D2_over_insulation_mm": 2e3 * self.R_B,
            "D3_outer_mm": 2e3 * self.R_C,
            "S_eff_mm2": self.effective_cross_section * 1e6,
            "S_geometric_mm2": self.conductor_geometric_area * 1e6,
            "fill_factor": self.effective_cross_section / self.conductor_geometric_area,
            "r20_ohm_per_km": float(self.resistance_per_length(20.0)) * 1e3,
            "conductor_mass_kg_per_m": self.conductor_mass,
            "conductor_density_effective": self.conductor_density_effective,
            "insulation": self.insulation_material.name,
            "jacket": self.jacket_material.name,
            "ambient_C": self.ambient_temperature,
        }


# ---------------------------------------------------------------------------
# Table 1 / Table 2 : Barrinax U-1000 R2V 1x70 mm² (XLPE / PVC)
# ---------------------------------------------------------------------------
CABLE_70 = CableSpec(
    name="cable70",
    designation="Barrinax U-1000 R2V 1x70 mm2 (0.6/1 kV, XLPE/PVC)",
    conductor_diameter=9.5e-3,
    insulation_thickness=1.1e-3,
    jacket_thickness=1.5e-3,
    effective_cross_section=70e-6,
    resistivity_T0=1.85e-8,
    alpha=0.0043,
    conductor_mass_per_length=0.584,
    ambient_temperature=19.0,
    test_steps=((120, 1950), (145, 1950), (170, 1500), (185, 1850), (212, 1350)),
    notes=(
        "All values from paper Table 1 and Table 2.",
        "XLPE and PVC densities are not given in the paper (assumed 920 / 1400 kg/m3).",
    ),
)

# ---------------------------------------------------------------------------
# Table 3 / Table 4 : H07RN-F TITANEX 1x150 mm²
# ---------------------------------------------------------------------------
CABLE_150 = CableSpec(
    name="cable150",
    designation="H07RN-F TITANEX 1x150 mm2",
    conductor_diameter=15e-3,
    insulation_thickness=2.0e-3,
    jacket_thickness=2.5e-3,
    effective_cross_section=150e-6,
    resistivity_T0=1.85e-8,
    alpha=0.0043,                 # not in Table 3 - same value as Table 1 assumed
    conductor_mass_per_length=1.74,
    ambient_temperature=23.5,
    test_steps=((130, 766), (163, 554), (220, 850), (305, 530), (405, 400)),
    notes=(
        "Values from paper Table 3 and Table 4.",
        "Table 3 gives neither the insulation/jacket materials nor alpha; the paper's XLPE/PVC "
        "property curves and alpha = 0.0043 1/K are used (a real H07RN-F has EPR insulation and a "
        "polychloroprene sheath).",
        "Table 3 conductor mass 1.74 kg/m exceeds the mass of 150 mm2 of solid copper (1.34 kg/m) "
        "and even of the full 15 mm conductor circle (1.58 kg/m); it is used as published.",
    ),
)

PRESETS = {"cable70": CABLE_70, "70": CABLE_70, "cable150": CABLE_150, "150": CABLE_150}


def get_cable(name: str) -> CableSpec:
    try:
        return PRESETS[str(name).lower()]
    except KeyError as exc:
        raise KeyError(f"unknown cable preset {name!r}; available: cable70, cable150") from exc
