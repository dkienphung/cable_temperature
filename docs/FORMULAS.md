# Formula check of Riba & Llauradó, *Materials* 2022, 15, 6814

"A Model to Calculate the Current–Temperature Relationship of Insulated and Jacketed Cables",
doi:10.3390/ma15196814 (CC BY 4.0).

Every equation of the paper (Eqs. 1–27) was transcribed from the article's MathML (Europe PMC
full text) and from the PDF, implemented without simplification in `cable_thermal/model.py`,
and checked by:

* **derivation / dimensional analysis** (each term in W/m or W/(m·K));
* **the cited source** (IEEE Std 738-2012 for Eqs. 25–26, CIGRE TB 207 / IEEE 738 for Eq. 1);
* **numerical tests** (`tests/test_model.py`, 17 tests): TDMA vs dense solver, exact mesh
  geometry, discrete energy conservation, steady state vs analytic cylinder solution, adiabatic
  short circuit vs IEC 60949, independent 121-node Radau reference solver, IEC 60287 consistency.

Legend: ✅ correct · ⚠️ correct but ambiguous/incomplete · ❌ error.

## 1. Governing equations

| Eq. | Content (as published) | Verdict | Notes |
|---|---|---|---|
| 1 | `I²_RMS r(T) = P_c + P_r − P_s + m c_p(T) dT/dt` [W/m] | ✅ | Heat balance of a **bare** conductor (CIGRE/IEEE 738). Moving the solar gain to the right gives −P_s – sign correct. Only used as background; implemented as `simulate_bare_conductor_eq1` for comparison (Fig. 17). |
| 2 | `r(T) = r_T0 [1 + α (T − T0)]` | ✅ / ⚠️ data | Formula correct. **But Table 1 gives α = 0.0043 K⁻¹ with T0 = 20 °C.** For annealed copper α₂₀ = 1/(234.5+20) = 0.00393 K⁻¹; 0.0043 ≈ 1/234.5 is the coefficient referred to **0 °C**. Effect on the 70 mm² test: +0.19 K (Fig. 16). The model uses the published 0.0043 (switchable). |
| 3–4 | `a_i T_{i−1}^{j+1} + b_i T_i^{j+1} + c_i T_{i+1}^{j+1} − d_i = 0`, tri-diagonal matrix with `a_1 = c_tot = 0` | ✅ | The printed matrix shows `a_1` / `c_tot` in the corners (cyclic layout); with `a_1 = c_tot = 0` it is the ordinary tri-diagonal system. Solved with the Thomas algorithm (test: identical to dense LU to 1e-13). |

## 2. Conductor (material A)

| Eq. | Content | Verdict | Notes |
|---|---|---|---|
| 5–7 | First node (centre): `ρ_A S_1,A c_p,A (T₁^{j+1}−T₁^j)/Δt = I_i² r_i(T₁^j) − k_A (T₁^{j+1}−T₂^{j+1}) p₁/Δx_A`, `S_1,A = π(x₁+0.5Δx_A)²`, `p₁ = 2π(x₁+0.5Δx_A)` | ✅ | x₁ = 0 → half element of radius Δx/2 (finite-volume control volume). For uniform generation this centre stencil is exact (T₀ − T(Δx) = q_v Δx²/4k). |
| 8–10 | Generic node: `S_i,A = π[(x_i+0.5Δx)² − (x_i−0.5Δx)²]`, `p_{i∓1} = 2π(x_i ∓ 0.5Δx)` | ✅ | `p_{i−1}`/`p_{i+1}` are the perimeters of the inner/outer faces of the control volume, so `k p/Δx` is the face conductance per unit length [W/(m·K)]. Coefficients (10) follow exactly. |
| 11–13 | Boundary node conductor/insulation, two half elements (A and B) | ✅ | `S_i,A = π[x_i² − (x_i−0.5Δx_A)²]`, `S_i,B = π[(x_i+0.5Δx_B)² − x_i²]` – checked: the half-areas of all nodes sum exactly to π R_A², π(R_B²−R_A²), π(R_C²−R_B²). Typo only: `r_i(T^j_{,i})`. |
| – | `J_i = I_i/S_i`, `r_i = ρ_e/S_i,A`, `I_i = S_i,A J_i` | ⚠️ | Correct (`I_i² r_i = J² ρ_e S_i,A`, i.e. volumetric Joule heat). The paper does not say which **total** area defines J: with J = I/(πR_A²) the total loss is I²ρ_e/(πR_A²) = I²·ρ_e/70.88 mm², with the "effective cross-section" (Table 1) it is I²ρ_e/70 mm². Implemented: homogenised conductor region so that Σ I_i² r_i = I² ρ_e/S_eff (Table 1); option `joule_area='geometric'` (−0.31 K). |
| – | Conductor heat capacity `ρ_A S_i,A c_p,A` | ⚠️ | Density of the stranded conductor region not given. Table 1 gives the conductor mass m = 0.584 kg/m, used as ρ_A = m/(πR_A²) so that Σρ_A S_i,A = m (consistent with Eq. 1). |

## 3. Insulation (B) and jacket (C)

| Eq. | Content | Verdict | Notes |
|---|---|---|---|
| 14, 16 | Generic insulation node | ✅ | |
| 15 | `S_i,A = π[(x_i+0.5Δx_B)² − (x_i−0.5Δx_B)²]` | ❌ typo | Must be **S_i,B** (it is the insulation area; Eq. 14 and 16 use S_i,B). Implemented as S_i,B. |
| 17–19 | Boundary node insulation/jacket | ✅ | |
| 20–22 | Generic jacket node | ✅ | |
| 23–24, 27 | Last node: `… = k_C (T_{i−1}−T_i)/Δx_C p_{i−1} − h p_{i+1}(T_i^{j+1} − T_air) − εσ p_{i+1}[(T_i^j)⁴ − T_air⁴]`, `p_{i+1} = 2πR_C` | ✅ | Radiation lagged (explicit) in time, convection implicit, temperatures in K for radiation. First-order semi-implicit; with Δt = 10 s the difference to fully implicit radiation is 0.01 K (Fig. 16). |
| – | k at a face | ⚠️ | Not specified for temperature-dependent k. Implemented with k at the mean temperature of the two nodes, which keeps `a_{i+1} = c_i` → the scheme is exactly conservative (discrete energy balance closes to 1e-9 W/m). |

## 4. Surface heat transfer

| Eq. | Content (as published) | Verdict | Correct form |
|---|---|---|---|
| 25 | `h = (3.645/π) ρ_air^0.25 D^−0.25 (T − T_air)^0.25` | ❌ | IEEE 738-2012 Eq. (5c): `q_cn = 3.645 ρ_f^0.5 D^0.75 (T_s−T_a)^1.25` [W/m]. Dividing by `πD(T_s−T_a)`: **`h = (3.645/π) ρ_f^0.5 D^−0.25 (T_s−T_a)^0.25`**. The exponent of ρ is 0.5, not 0.25. With ρ ≈ 1.18 kg/m³ the printed formula gives h about 4 % too low (70 mm² test: +0.42 K). Both versions are available (`--convection paper|ieee738`); the corrected one is the default. |
| 26 | `ρ_air = (1.293 − 1.525·10⁻⁴H + 6.379·10⁻⁹H²) / (1 + 0.00367 (T_air+T)/2)` with "T_air (K) and T (K)" | ❌ units | The formula is IEEE 738-2012 Eq. (14a) and requires the **film temperature in °C** (0.00367 ≈ 1/272.5). Taken literally in kelvin it gives ρ ≈ 0.62 kg/m³ at 21 °C (instead of 1.20), i.e. h −28 %. Implemented in °C. The elevation H is not given in the paper (H = 0 used; 280 m → +0.18 K). |
| – | ε = 0.85, σ | ✅ | σ = 5.670374419·10⁻⁸ W/(m²K⁴). |

## 5. Data (Tables 1–4) and internal consistency

| Item | Finding |
|---|---|
| Table 1, ρ = 1.85·10⁻⁸ Ω·m, S = 70 mm² | r₂₀ = 0.264 Ω/km – plausible (IEC 60228 class 2 limit 0.268 Ω/km). |
| Table 1, m = 0.584 kg/m | = 65.2 mm² of copper – plausible for a 70 mm² stranded conductor. |
| Table 3, m = 1.74 kg/m (150 mm²) | ❌ physically impossible for the conductor: 150 mm² of copper weigh 1.34 kg/m and even the full 15 mm circle of solid copper only 1.58 kg/m. Probably the mass of the complete cable. Used as published (flagged). |
| Table 3 materials / α | Not given. H07RN-F is a harmonised 450/750 V cable with **rubber (EPR) insulation and polychloroprene sheath**, not XLPE/PVC; the paper's XLPE/PVC curves and α = 0.0043 are used, as nothing else is given. |
| Densities of XLPE / PVC | Not given (needed in Eqs. 13–27). Assumed 920 and 1400 kg/m³; ±20 % changes the end temperature by ∓0.15 K. |
| Property curves (Figs. 2–3) | Only cover 20–75 °C. Digitized (tick-calibrated, < 0.3 % error) and held constant outside – relevant for 90 °C ratings and short circuits. |
| Fig. 10b Joule heating | ❌ inconsistent with Table 1: I²r(T) from Table 1 data is 3.79 W/m at 120 A/19 °C and 12.8 W/m at 212 A/≈37 °C; Fig. 10b shows ≈3.4 and ≈10.7 W/m (−10 % … −16 %). |
| Fig. 10b losses vs Fig. 9 temperatures | ❌ inconsistent with Eqs. 23–26: at the Fig. 9 surface temperatures, P_c + P_r from Eqs. 25–26 + radiation are 2.0 W/m (t = 1950 s) and 9.5 W/m (t = 8550 s); Fig. 10b shows 3.0 and 10.7 W/m. |
| Fig. 10a | 28 nodal elements = 9+9+9+1 nodes – consistent with "nine nodes for each material". |

## 6. Consequence: can Figs. 9–11 be reproduced?

* **Numerics are not the issue.** The paper's scheme (9 nodes per material, Δt = 10 s) differs
  from an independent 121-node Radau solution by < 0.01 K for both tests.
* **150 mm² cable (Fig. 11):** reproduced with the published data – mean error −0.4 K at the
  conductor-insulation and insulation-jacket points (RMS 0.46 / 0.40 K, mean relative error
  −1.5 %), jacket-air −1.0 K.
* **70 mm² cable (Fig. 9):** **not reproduced** with the published data – the model is
  +3.2 K (mean) / +5.3 K (end) hotter. No formula correction explains this (all corrections
  together < 1 K, Fig. 16). The early heating rate in Fig. 9 (≈ 2.5 K in the first 500 s with
  ≈ 3.5–3.8 W/m) requires about 2.3× the heat capacity per unit length obtainable from Table 1
  (368 J/(m·K)), and the late level about 1.3× the surface heat-transfer coefficient
  (diagnostic least-squares fit, `results/RESULTS.md` §2). Together with the Fig. 10b
  inconsistencies above, this indicates that the published Fig. 9/10 simulation used inputs
  different from Table 1 (or a different heat capacity/heat-transfer implementation than
  printed).

## 7. Independent checks of the implementation

| Check | Result |
|---|---|
| Steady state, constant k, vs exact cylinder solution `ΔT = q/(2πk) ln(r₂/r₁)`, `q/(4πk_A)` | < 2 % of ΔT with 9 nodes, < 0.2 % with 40 nodes (2nd-order convergence) |
| Discrete energy balance Σ C_i ΔT_i/Δt = ΣI²r − P_c − P_r | 5·10⁻¹⁰ W/m |
| Adiabatic short circuit vs IEC 60949 `I²t = K²S² ln((θf+β)/(θi+β))` | < 0.05 K (K = 226 recovered from σ_c = 3.45 MJ/(K·m³), β = 234.5 K) |
| Non-adiabatic short circuit (70 mm², 9.3 kA, 1 s, from 90 °C) | model 233.8 °C vs IEC 60949 non-adiabatic factor 231.9 °C |
| Steady-state ampacity at 90 °C | model 339 A vs IEC 60287 (same k) 354 A (70 mm²); 550 A vs 561 A (150 mm²) – IEC uses its empirical h = Z/De^g + E |
| Long transient → steady state | identical to 1e-3 K |
