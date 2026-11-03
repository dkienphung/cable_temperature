# Results - cable current-temperature model (Riba & Llauradó, Materials 2022, 15, 6814)

All numbers below are produced by `python -m cable_thermal reproduce` (file generated automatically).

Paper curves were digitized from the published figures (±0.15 °C reading accuracy). For the conductor-insulation and insulation-jacket traces the paper's measured and simulated curves overlap and the digitized trace is their median; the jacket-air trace is the measured one.

## 1. Test sequences of the paper (Figs. 9 and 11)

Model = paper Eqs. (2)-(27) with the cable data of Tables 1-4, 9 nodes per material, Δt = 10 s, site elevation H = 280 m (Terrassa laboratory, Eq. 26), still air as assumed by the paper.

### 70 mm² (Table 1/2, Fig. 9)

| Measuring point | paper end [°C] | model end, Eq. 25 corrected [°C] | model end, Eq. 25 as printed [°C] | mean error [K] | RMS error [K] | mean rel. error [%] |
|---|---|---|---|---|---|---|
| Conductor-insulation | 39.31 | 44.77 | 45.10 | +3.26 | 3.40 | +10.6 |
| Insulation-jacket | 37.52 | 43.47 | 43.80 | +3.80 | 3.96 | +12.8 |
| Jacket-air | 35.21 | 40.12 | 40.44 | +3.05 | 3.16 | +10.8 |

### 150 mm² (Table 3/4, Fig. 11)

| Measuring point | paper end [°C] | model end, Eq. 25 corrected [°C] | model end, Eq. 25 as printed [°C] | mean error [K] | RMS error [K] | mean rel. error [%] |
|---|---|---|---|---|---|---|
| Conductor-insulation | 39.54 | 38.68 | 38.74 | -0.42 | 0.45 | -1.4 |
| Insulation-jacket | 37.85 | 37.25 | 37.31 | -0.36 | 0.39 | -1.2 |
| Jacket-air | 37.00 | 34.50 | 34.58 | -0.95 | 1.13 | -3.2 |

Discrete energy balance error of the TDMA solution (max over all steps, W/m): cable70/paper: 7.0e-10, cable70/ieee738: 6.2e-10, cable150/paper: 5.4e-10, cable150/ieee738: 4.6e-10

Heat balance at t = 8600 s (70 mm², model): I²r = 13.14 W/m, P_c = 7.40 W/m, P_r = 5.22 W/m, surface 40.12 °C. Paper Fig. 10b: I²r ≈ 10.8, P_c ≈ 6.3, P_r ≈ 4.4 W/m.

## 2. Diagnostic inverse analysis

Least-squares multipliers on the total heat capacity and on the convection coefficient that make the model reproduce the digitized paper curves (diagnosis only - not used elsewhere):

| Cable | heat-capacity factor | h factor | RMS after fit [K] |
|---|---|---|---|
| cable70 | 2.34 | 1.35 | 0.41 |
| cable150 | 0.85 | 1 (not fitted) | 0.41 |

## 3. Numerical verification

* cable70: paper scheme vs independent 121-node Radau reference, max |ΔT| = 0.0068 K. Convergence: N=3 (Δt=40 s): 0.1183 K, N=5 (Δt=20 s): 0.0912 K, N=9 (Δt=10 s): 0.0068 K, N=18 (Δt=5 s): 0.0036 K, N=36 (Δt=2.5 s): 0.0019 K.
* cable150: paper scheme vs independent 121-node Radau reference, max |ΔT| = 0.0086 K. Convergence: N=3 (Δt=40 s): 0.2343 K, N=5 (Δt=20 s): 0.0727 K, N=9 (Δt=10 s): 0.0086 K, N=18 (Δt=5 s): 0.0099 K, N=36 (Δt=2.5 s): 0.0100 K.

## 4. Steady-state current-temperature relationship (ampacity at 90 °C)

| Cable | T_air [°C] | model [A] | IEC 60287, k from paper curves [A] | IEC 60287, IEC ρ_T [A] |
|---|---|---|---|---|
| cable70 | 19 | 338 | 354 | 372 |
| cable150 | 23.5 | 548 | 561 | 596 |

The paper's property curves end at 75 °C and are held constant above (documented limitation).

## 5. Short circuit (70 mm²)

Load 338 A (steady 90.0 °C), then 9.27 kA for 1 s (the current that IEC 60949 predicts gives 250 °C adiabatically with the same copper data, K = 203.4 A·s^0.5/mm², β = 212.6 K; standard IEC copper K = 225.7).

* model peak conductor temperature: 233.8 °C (centre), 232.8 °C (conductor-insulation boundary)
* IEC 60949 adiabatic: 250.0 °C (upper bound: no heat leaves the conductor)
* IEC 60949 non-adiabatic (ε = 1 + X√(t/S) + Y t/S, Cu/XLPE X = 0.41, Y = 0.12 → ε = 1.051): 231.9 °C - the model's heat flow into the XLPE gives practically the same result
* with the paper's Δt = 10 s the 1 s pulse is not resolved (peak 90.0 °C) - short circuits need Δt ≪ t_sc.

## 6. Sensitivity of the 70 mm² end-of-test temperature

Base (Table 1 data, Eq. 25 corrected): 44.77 °C at t = 8600 s; paper: 39.31 °C.

| Variant | ΔT [K] | RMS vs paper [K] |
|---|---|---|
| Eq. (25) as printed (ρ^0.25) | +0.33 | 3.72 |
| α = 0.00393 1/K (Cu @ 20 °C) | -0.19 | 3.45 |
| ρ_e = 1.7241e-8 Ω·m (IEC 60028) | -1.70 | 2.61 |
| J over π R_A² (joule_area=geometric) | -0.31 | 3.36 |
| ε = 0.95 | -0.86 | 3.02 |
| ε = 0.75 | +0.94 | 4.07 |
| Elevation H = 0 m (sea level) | -0.18 | 3.43 |
| Indoor air speed 0.1 m/s | -1.68 | 2.05 |
| XLPE & PVC density +20 % | -0.15 | 3.39 |
| XLPE & PVC density −20 % | +0.15 | 3.66 |
| Conductor mass = 8960·S_eff | -0.09 | 3.44 |
| Constant properties @20 °C | -0.02 | 3.55 |
| Radiation linearized (implicit) | -0.01 | 3.51 |
| Fig. 10b I²r level (joule ×0.90) | -2.50 | 2.19 |

## 7. Dynamic case and Eq. (1)

Peak conductor-insulation temperature over 48 h: 48.7 °C (variable ambient) vs 42.9 °C (constant 19 °C).

Bare conductor with lumped Eq. (1): 48.87 °C at 8600 s.

## 8. Laboratory environment: elevation and indoor air movement

The paper does not state the elevation (Eq. 26) and assumes still air. The paper cables are preset to the Terrassa laboratory (H = 280 m). Optional indoor air movement uses IEEE 738-2012 forced convection (largest of natural, q_c1, q_c2). Conductor-insulation temperature at the end of each test:

| Case | 70 mm² end [°C] | 70 mm² RMS [K] | 150 mm² end [°C] | 150 mm² RMS [K] |
|---|---|---|---|---|
| Paper | 39.31 | - | 39.54 | - |
| H = 0 m, still air | 44.59 | 3.43 | 38.65 | 0.74 |
| H = 280 m, still air | 44.77 | 3.52 | 38.68 | 0.73 |
| H = 280 m, V = 0.05 m/s | 44.77 | 3.47 | 38.62 | 0.79 |
| H = 280 m, V = 0.1 m/s | 43.09 | 2.05 | 38.05 | 1.06 |
| H = 280 m, V = 0.15 m/s | 40.80 | 0.94 | 37.58 | 1.27 |
| H = 280 m, V = 0.2 m/s | 39.21 | 0.79 | 37.23 | 1.43 |
| H = 280 m, V = 0.25 m/s | 38.03 | 1.23 | 36.94 | 1.56 |
| H = 280 m, V = 0.3 m/s | 37.09 | 1.70 | 36.70 | 1.67 |
| H = 280 m, V = 0.4 m/s | 35.68 | 2.46 | 36.32 | 1.85 |
| H = 280 m, V = 0.5 m/s | 34.65 | 3.04 | 36.01 | 1.99 |

Best-fitting air speed: 70 mm² 0.18 m/s (RMS 0.73 K), 150 mm² 0.00 m/s (RMS 0.73 K). The elevation changes the results by only ~0.2 K; a gentle draught of ~0.2 m/s reconciles the 70 mm² test but makes the 150 mm² test worse, so a single air speed does not explain both tests.
