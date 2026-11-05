# Cable temperature model

A Python application that implements the transient thermal model of insulated and jacketed
cables from

> J.-R. Riba, J. Llauradó, *A Model to Calculate the Current–Temperature Relationship of Insulated
> and Jacketed Cables*, **Materials** 2022, 15(19), 6814. https://doi.org/10.3390/ma15196814

The model follows the paper exactly, with nothing simplified (Eqs. 1–27):

* 1-D radial finite differences: copper conductor (A), XLPE insulation (B), PVC jacket (C), with
  interface nodes. Nine nodes per material gives 28 nodes, as in the paper.
* Implicit time stepping (Δt = 10 s) of the TDMA system a·T₍i−1₎ + b·T₍i₎ + c·T₍i+1₎ = d. Every
  coefficient is coded exactly as in Eqs. 5–27, and the system is solved with the Thomas algorithm.
* Temperature-dependent resistance r(T) (Eq. 2). Temperature-dependent c_p(T) and k(T) of Cu,
  XLPE and PVC, digitized from the paper's Figs. 2–3.
* Natural convection (IEEE 738, Eqs. 25–26) and radiation with ε = 0.85 (Eq. 23).

On top of the paper, every formula was checked and the solver was verified independently. The
issues found are in **[docs/FORMULAS.md](docs/FORMULAS.md)**; the numerical results are in
**[results/RESULTS.md](results/RESULTS.md)**.

The model is also applied to a **database of 337 real manufacturer wires (36 AWG to 1/0 AWG;
PVC, silicone rubber, PTFE, XLPE, mPPE, ETFE, PVDF, EPR, fluoroplastic)** taken from Alpha Wire
datasheets, and compared with manufacturer current ratings (Habia, RS PRO, Alpha Wire, General
Cable): **[docs/WIRE_DATABASE.md](docs/WIRE_DATABASE.md)**,
**[results/wires/RESULTS_WIRES.md](results/wires/RESULTS_WIRES.md)**.

**MATLAB / Octave / Simulink version:** the folder [`matlab/`](matlab/README.md) contains a
self-contained port: the model, the paper study, the wire database study, tests with Python
parity, and a Simulink model builder.

## Quick start

```bash
pip install -r requirements.txt
```

```bash
python -m cable_thermal reproduce --out results
```

```bash
python -m cable_thermal simulate --cable cable70 --plot run.png --csv run.csv
```

```bash
python -m cable_thermal simulate --cable cable70 --steps 150:3600,250:1800 --ambient 25
```

```bash
python -m cable_thermal simulate --config examples/overload.json
```

```bash
python -m cable_thermal steady --cable cable150 --current 400 --ambient 30
```

```bash
python -m cable_thermal ampacity --cable cable70 --tmax 90 --ambient 30
```

```bash
python -m pytest -q
```

Options:
* `--nodes`: nodes per material.
* `--dt`: time step.
* `--convection ieee738|paper`: corrected vs printed Eq. 25.
* `--elevation`: H in Eq. 26. The default is the cable preset; the paper's two cables are preset
  to 280 m, the Terrassa lab.
* `--air-speed` / `--air-angle`: indoor air movement, see below.

A JSON scenario file can set the same options under `"model"`: `elevation`, `air_speed`,
`air_angle`, `nodes`, `dt` and `convection`.

From Python:

```python
from cable_thermal import CABLE_70, CableThermalModel
m = CableThermalModel(CABLE_70)                      # Table 1 data, 9 nodes/material, dt = 10 s
r = m.simulate(CABLE_70.test_steps)                  # Table 2 current steps
r.T_conductor_insulation, r.T_insulation_jacket, r.T_jacket_air, r.joule, r.p_conv, r.p_rad
m.steady_state(212.0); m.ampacity(90.0)
CableThermalModel(CABLE_70, air_speed=0.2, elevation=0.0)   # draught of 0.2 m/s, sea level
m.simulate(lambda t: 9270 if t < 1 else 339, t_end=60, dt=0.005, T_init=m.steady_state(339))  # short circuit
```

## Formula check: summary (details in docs/FORMULAS.md)

| Item | Finding |
|---|---|
| Eqs. 1–24, 27 (heat balance, FDM/TDMA coefficients, boundary nodes) | Correct. Areas, perimeters and coefficients check out, the scheme conserves energy (residual 5·10⁻¹⁰ W/m), and it is 2nd-order convergent. |
| **Eq. 25** natural convection | ❌ `ρ_air^0.25` must be `ρ_air^0.5` (IEEE 738-2012: q_cn = 3.645 ρ^0.5 D^0.75 ΔT^1.25). |
| **Eq. 26** air density | ❌ the film temperature must be in °C, not K as the text says (otherwise ρ comes out about 2× too low). |
| **Eq. 15** | ❌ typo: S_i,A should be S_i,B. |
| Table 1: α = 0.0043 K⁻¹ at T0 = 20 °C | ⚠️ That is the 0 °C value; copper at 20 °C has α = 0.00393 K⁻¹. |
| Table 3: conductor mass 1.74 kg/m (150 mm²) | ❌ More than 150 mm² of solid copper weighs (1.34 kg/m). H07RN-F is EPR/PCP, not XLPE/PVC. |
| Not given | XLPE/PVC densities, elevation H, how J and k are evaluated at faces. Assumptions are documented. |
| Fig. 10b | ❌ I²r is 10–16 % below Table 1 data, and P_c + P_r does not match Eqs. 23–26 at the Fig. 9 temperatures. |

## Results (figures in `results/figures`)

| Figure | Content |
|---|---|
| `fig02_specific_heat.png`, `fig03_thermal_conductivity.png` | Digitized property curves (paper Figs. 2–3) |
| `fig04_discretization.png` | Radial mesh (paper Figs. 4–7) |
| `fig09_cable70_temperature_rise.png` | 70 mm², 5 current steps: model vs paper Fig. 9, plus relative difference |
| `fig09b_cable70_formula_variants.png` | Eq. 25 printed vs corrected, plus the diagnostic fit |
| `fig10a_cable70_radial_profiles.png` | Radial temperature profile every 1000 s (paper Fig. 10a) |
| `fig10b_cable70_heat_balance.png` | Joule, P_c, P_r, P_c + P_r (paper Fig. 10b) |
| `fig11_cable150_temperature_rise.png` | 150 mm², 5 current steps: model vs paper Fig. 11 |
| `fig12_verification_reference_solver.png` | Paper scheme vs independent Radau solver; convergence |
| `fig13_steady_state_current_temperature.png` | Steady-state current–temperature curves, ampacity at 90 °C, IEC 60287 |
| `fig14_short_circuit.png` | 1 s short circuit vs IEC 60949 (adiabatic and non-adiabatic) |
| `fig15_dynamic_current_ambient.png` | 48 h with a daily load cycle, an overload and variable ambient temperature |
| `fig16_sensitivity.png` | Sensitivity to formula choices and missing data |
| `fig17_bare_conductor_eq1.png` | Insulated cable vs bare-conductor Eq. 1 |
| `fig18_heat_transfer_coefficients.png` | h from Eq. 25 printed/corrected, Eq. 26 in K, radiation |
| `fig19_lab_environment_air_speed.png` | Lab elevation and indoor air movement vs the paper curves |

Main numbers:

All results use the Terrassa lab elevation (280 m) and still air, as in the paper.

* **150 mm² (Fig. 11)** is reproduced from the published data. At the end of the test the model
  gives 38.7 / 37.3 / 34.5 °C against the paper's 39.5 / 37.9 / 37.0 °C. Mean error is −0.4 K
  at the conductor.
* **70 mm² (Fig. 9)** is not reproduced from the published data. The model ends at
  44.8 / 43.5 / 40.1 °C against the paper's 39.3 / 37.5 / 35.2 °C. None of the formula
  corrections explains the gap, and the numerics are exact to 0.01 K. An inverse fit needs about
  2.3× the heat capacity in Table 1 and 1.3× the heat-transfer coefficient, which points to
  different inputs behind the published simulation (see FORMULAS.md §6).
* Over the 20–75 °C range of the tests, using temperature-dependent properties instead of
  constant ones (at 20 °C) changes the result by only −0.02 K.
* Ampacity at 90 °C: 338 A for 70 mm² at 19 °C (IEC 60287 with the same k: 354 A) and 548 A
  for 150 mm² at 23.5 °C (IEC 60287: 561 A).
* Short circuit, 70 mm², 9.27 kA for 1 s from 90 °C: peak 233.8 °C. IEC 60949 gives 250 °C
  adiabatic and 231.9 °C non-adiabatic. The paper's Δt = 10 s cannot resolve this case.

## Laboratory environment: elevation and indoor air movement

The paper does not state the elevation of its lab and assumes still air. Two options cover this:

* **Site elevation** (`CableSpec.elevation`, `--elevation`) enters the air density of Eq. 26.
  The paper's cables are preset to 280 m (UPC Terrassa). Compared with sea level this changes
  the results by only about +0.2 K, because thinner air cools slightly less.
* **Indoor air speed** (`ModelOptions.air_speed`, `--air-speed`, `--air-angle`) is not part of the
  paper. It adds IEEE 738-2012 forced convection (Eqs. 3a/3b, with the air properties of
  Eqs. 13a/15a and the direction factor of Eq. 4a). As IEEE 738 prescribes, the larger of the
  natural and forced losses is used, so `air_speed = 0` reproduces the paper exactly.

Results (`fig19`, RESULTS.md §8): a gentle draught of about 0.18 m/s brings the 70 mm² test from
3.5 K RMS down to 0.7 K RMS against the paper (end 39.2 °C vs 39.3 °C). The same draught makes
the 150 mm² test worse, and that test fits best at 0 m/s, so a single air speed cannot explain
both tests. The too-fast early heating of the 70 mm² cable, which is driven by heat capacity, is
also not affected by air movement.

## Wire database (36 AWG to 1/0 AWG)

```bash
python -m cable_thermal wire-list --awg 22
```

```bash
python -m cable_thermal wire --part 391250 --ambient 30 --current 8 --duration 600 --plot w.png
```

```bash
python -m cable_thermal wire-study --out results/wires
```

* **Data:** 337 Alpha Wire hook-up wire parts parsed from the manufacturer's spec sheets: AWG,
  stranding, plating, conductor and outer diameter, insulation material and wall, jacket,
  temperature rating, voltage and DC resistance. Insulations covered: PVC 123, silicone 73,
  PTFE 58, mPPE 24, XLPE 20, ETFE 11, fluoroplastic 11, PVDF 4, EPR 3, and mica-glass 10 (the
  mica-glass parts are not simulated, as no thermal data exist).
* **Reference ratings:** Habia PTFE type E (25 types), RS PRO PTFE type A (3), Alpha PV wire (3)
  and Super Vu-Tron welding cable (6 AWG to 4/0, 8 sizes).
* **Material properties:** k, ρ and c_p for each insulation, with the source of every value
  (`cable_thermal/data/insulation_materials.json`).
* **Convection:** for thin wires the paper's Eq. 25 is outside its range of validity, so the wire
  study uses Churchill & Chu (option `--convection churchill_chu`). Eq. 25 under-predicts the
  ampacity by about 16 % for wires under 1 mm OD.

Comparison of the model with the manufacturers' ratings:

| Data set | Model result |
|---|---|
| RS PRO PTFE, rated for a 40 K rise | 35–43 K rise at the published currents |
| Alpha PV wire and Super Vu-Tron welding cable (90 °C, 8 AWG to 4/0) | 76–78 °C at the published currents, so the ratings are 0.87–0.89 × the model's 90 °C ampacity |
| Habia PTFE (200 °C) | conductor at 112–191 °C at the published currents, so the ratings stay below the 200 °C limit |

Figures are `results/wires/figures/fig20`–`fig26`; per-part results are in
`results/wires/wire_results.csv`.

## Repository layout

```
cable_thermal/      model.py (paper FDM/TDMA), cables.py (Tables 1-4), properties.py,
                    heat_transfer.py (Eqs. 23-26 + forced convection + Churchill-Chu), iec.py
                    (IEC 60287 / 60949), reference.py (independent solver), reproduce.py (all
                    figures), wires.py (wire database -> model), wire_study.py (wire campaign),
                    __main__.py (CLI)
cable_thermal/data  digitized property curves, paper result curves, insulation materials (JSON)
cable_thermal/data/wires  Alpha Wire database (337 parts) and reference current ratings (JSON)
tools/              digitizers for the paper's figures; Alpha Wire fetch/parse; reference parsers
docs/               FORMULAS.md (equation check), WIRE_DATABASE.md (wire data, sources, method)
matlab/             MATLAB / Octave port (core, heat, iec, wires, simulink, run_*.m scripts)
results/            figures, RESULTS.md, summary.json, CSV time series
results/wires/      wire campaign: RESULTS_WIRES.md, wire_results.csv, comparisons.json, figures
tests/              verification tests
```

The paper's curves were digitized from its published figures. The tick checks agree to better
than 0.3 % of the axis span, and readings are accurate to about ±0.15 °C. For the
conductor–insulation and insulation–jacket points, the paper's measured and simulated curves lie
on top of each other, so the digitized curve is the median of the two. The paper is CC BY 4.0.
