# Results - cable current-temperature model (Riba & Llauradó, Materials 2022, 15, 6814)

All numbers below are produced by `python -m cable_thermal reproduce` (file generated automatically).

Paper curves were digitized from the published figures (±0.15 °C reading accuracy). For the conductor-insulation and insulation-jacket traces the paper's measured and simulated curves overlap and the digitized trace is their median; the jacket-air trace is the measured one.

## 1. Test sequences of the paper (Figs. 9 and 11)

Model = paper Eqs. (2)-(27) with the cable data of Tables 1-4, 9 nodes per material, Δt = 10 s.

### 70 mm² (Table 1/2, Fig. 9)

| Measuring point | paper end [°C] | model end, Eq. 25 corrected [°C] | model end, Eq. 25 as printed [°C] | mean error [K] | RMS error [K] | mean rel. error [%] |
|---|---|---|---|---|---|---|
| Conductor-insulation | 39.31 | 44.59 | 45.01 | +3.17 | 3.31 | +10.4 |
| Insulation-jacket | 37.52 | 43.29 | 43.71 | +3.71 | 3.86 | +12.6 |
| Jacket-air | 35.21 | 39.94 | 40.35 | +2.96 | 3.06 | +10.5 |

### 150 mm² (Table 3/4, Fig. 11)

| Measuring point | paper end [°C] | model end, Eq. 25 corrected [°C] | model end, Eq. 25 as printed [°C] | mean error [K] | RMS error [K] | mean rel. error [%] |
|---|---|---|---|---|---|---|
| Conductor-insulation | 39.54 | 38.65 | 38.73 | -0.43 | 0.46 | -1.5 |
| Insulation-jacket | 37.85 | 37.21 | 37.29 | -0.37 | 0.40 | -1.3 |
| Jacket-air | 37.00 | 34.46 | 34.56 | -0.96 | 1.15 | -3.2 |

Discrete energy balance error of the TDMA solution (max over all steps, W/m): cable70/paper: 6.9e-10, cable70/ieee738: 5.5e-10, cable150/paper: 4.1e-10, cable150/ieee738: 5.1e-10

Heat balance at t = 8600 s (70 mm², model): I²r = 13.13 W/m, P_c = 7.45 W/m, P_r = 5.17 W/m, surface 39.94 °C. Paper Fig. 10b: I²r ≈ 10.8, P_c ≈ 6.3, P_r ≈ 4.4 W/m.

## 2. Diagnostic inverse analysis

Least-squares multipliers on the total heat capacity and on the convection coefficient that make the model reproduce the digitized paper curves (diagnosis only - not used elsewhere):

| Cable | heat-capacity factor | h factor | RMS after fit [K] |
|---|---|---|---|
| cable70 | 2.34 | 1.33 | 0.41 |
| cable150 | 0.85 | 1 (not fitted) | 0.41 |

## 3. Numerical verification

* cable70: paper scheme vs independent 121-node Radau reference, max |ΔT| = 0.0068 K. Convergence: N=3 (Δt=40 s): 0.1183 K, N=5 (Δt=20 s): 0.0912 K, N=9 (Δt=10 s): 0.0068 K, N=18 (Δt=5 s): 0.0036 K, N=36 (Δt=2.5 s): 0.0019 K.
* cable150: paper scheme vs independent 121-node Radau reference, max |ΔT| = 0.0087 K. Convergence: N=3 (Δt=40 s): 0.2337 K, N=5 (Δt=20 s): 0.0725 K, N=9 (Δt=10 s): 0.0087 K, N=18 (Δt=5 s): 0.0099 K, N=36 (Δt=2.5 s): 0.0100 K.

## 4. Steady-state current-temperature relationship (ampacity at 90 °C)

| Cable | T_air [°C] | model [A] | IEC 60287, k from paper curves [A] | IEC 60287, IEC ρ_T [A] |
|---|---|---|---|---|
| cable70 | 19 | 339 | 354 | 372 |
| cable150 | 23.5 | 550 | 561 | 596 |

The paper's property curves end at 75 °C and are held constant above (documented limitation).

## 5. Short circuit (70 mm²)

Load 339 A (steady 90.0 °C), then 9.27 kA for 1 s (the current that IEC 60949 predicts gives 250 °C adiabatically with the same copper data, K = 203.4 A·s^0.5/mm², β = 212.6 K; standard IEC copper K = 225.7).

* model peak conductor temperature: 233.8 °C (centre), 232.8 °C (conductor-insulation boundary)
* IEC 60949 adiabatic: 250.0 °C (upper bound: no heat leaves the conductor)
* IEC 60949 non-adiabatic (ε = 1 + X√(t/S) + Y t/S, Cu/XLPE X = 0.41, Y = 0.12 → ε = 1.051): 231.9 °C - the model's heat flow into the XLPE gives practically the same result
* with the paper's Δt = 10 s the 1 s pulse is not resolved (peak 90.0 °C) - short circuits need Δt ≪ t_sc.

## 6. Sensitivity of the 70 mm² end-of-test temperature

Base (Table 1 data, Eq. 25 corrected): 44.59 °C at t = 8600 s; paper: 39.31 °C.

| Variant | ΔT [K] | RMS vs paper [K] |
|---|---|---|
| Eq. (25) as printed (ρ^0.25) | +0.42 | 3.67 |
| α = 0.00393 1/K (Cu @ 20 °C) | -0.19 | 3.36 |
| ρ_e = 1.7241e-8 Ω·m (IEC 60028) | -1.69 | 2.52 |
| J over π R_A² (joule_area=geometric) | -0.31 | 3.26 |
| ε = 0.95 | -0.85 | 2.94 |
| ε = 0.75 | +0.92 | 3.96 |
| Elevation H = 280 m (Terrassa) | +0.18 | 3.52 |
| XLPE & PVC density +20 % | -0.15 | 3.29 |
| XLPE & PVC density −20 % | +0.14 | 3.56 |
| Conductor mass = 8960·S_eff | -0.09 | 3.34 |
| Constant properties @20 °C | -0.02 | 3.45 |
| Radiation linearized (implicit) | -0.01 | 3.42 |
| Fig. 10b I²r level (joule ×0.90) | -2.48 | 2.10 |

## 7. Dynamic case and Eq. (1)

Peak conductor-insulation temperature over 48 h: 48.5 °C (variable ambient) vs 42.7 °C (constant 19 °C).

Bare conductor with lumped Eq. (1): 48.59 °C at 8600 s.
