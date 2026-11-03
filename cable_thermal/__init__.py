"""Cable current-temperature model after Riba & Llauradó, Materials 2022, 15, 6814."""
from .cables import CABLE_70, CABLE_150, CableSpec, get_cable
from .model import CableThermalModel, ModelOptions, SimulationResult, StepProfile, tdma

__all__ = ["CABLE_70", "CABLE_150", "CableSpec", "get_cable", "CableThermalModel", "ModelOptions",
           "SimulationResult", "StepProfile", "tdma"]
__version__ = "1.0.0"
