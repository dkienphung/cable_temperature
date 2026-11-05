# MATLAB / Octave / Simulink version of the cable temperature study

This folder is a self-contained MATLAB port of the Python package (`cable_thermal`). It covers:

* the model of Riba & Llauradó, *Materials* 2022, 15, 6814, implemented exactly as published:
  radial FDM and TDMA solution of Eqs. 2–27, temperature-dependent material properties
  digitized from the paper's Figs. 2–3, and convection and radiation per Eqs. 23–26;
* the formula checks: Eq. 25 printed vs corrected, Eq. 26 units, IEC 60287 and IEC 60949;
* the optional extensions: site elevation, indoor air speed (IEEE 738 forced convection) and
  Churchill–Chu natural convection for thin wires;
* the wire database: 337 Alpha Wire parts from 36 AWG to 1/0 AWG, plus the Habia, RS PRO,
  Alpha PV wire and Super Vu-Tron current ratings;
* a builder for a Simulink model.

Requirements: MATLAB R2016b or newer (for `jsondecode` and `strjoin`), with no toolboxes.
GNU Octave 7 or newer also works, and all scripts were run with Octave 8.4. For headless
Octave, run the scripts under `xvfb-run octave --no-gui` so that the Qt renderer is used; without a
display the code falls back to gnuplot, which renders the figures less well.

## Quick start

```matlab
cd matlab
setup_paths                 % add the folders to the path
run_tests                   % verification + parity with the Python results (all PASS)
run_paper_study             % paper study -> results/paper/
run_wire_study              % whole wire database -> results/wires/
```

Single runs:

```matlab
c = ct_cable_preset('cable70');                 % Table 1/2 (Barrinax U-1000 R2V 1x70 mm2)
m = ct_model(c);                                % 9 nodes per material, dt = 10 s, Eq. 25 corrected
r = ct_simulate(m, c.test_steps);               % current steps of Table 2
plot(r.t, ct_probe(r, 'conductor_insulation'))
T = ct_steady_state(m, 212);                    % steady radial temperature at 212 A
I = ct_ampacity(m, 90, 30);                     % current for 90 degC in 30 degC air

m2 = ct_model(c, [], 'convection', 'paper');    % Eq. 25 exactly as printed (rho^0.25)
m3 = ct_model(c, [], 'air_speed', 0.2);         % indoor draught 0.2 m/s (IEEE 738)

w  = ct_wire_find('3051');                      % Alpha Wire 22 AWG UL1007 PVC
mw = ct_model(ct_wire_to_cable(w, 30), [], 'convection', 'churchill_chu');
Iw = ct_ampacity(mw, w.temp_max_C, 30);         % 11.61 A
```

## Folder layout

| Folder | Contents |
|---|---|
| `core/` | model: `ct_model` (mesh of Figs. 4–7, areas and perimeters of Eqs. 6–24), `ct_coefficients` (TDMA coefficients of Eqs. 7–27), `ct_tdma` (Thomas algorithm), `ct_simulate`, `ct_steady_state`, `ct_ampacity`, `ct_time_constant`, materials, cable presets, options |
| `heat/` | Eq. 25 printed/corrected, Eq. 26 (air density), Churchill–Chu, IEEE 738 forced convection, radiation |
| `iec/` | IEC 60287 free-air rating, IEC 60949 adiabatic short circuit |
| `wires/` | AWG/stranding utilities, Alpha Wire / Habia / RS / Super Vu-Tron record conversion, wire → cable |
| `util/` | JSON loading, plotting/export helpers (MATLAB and Octave), paper-curve comparison |
| `data/` | the same JSON data as the Python package (digitized paper curves, materials, wire database, reference ratings) |
| `reference/` | Python results used by `run_tests` for the parity checks |
| `simulink/` | `build_cable_simulink_model.m`, `ct_simulink_step.m`, `run_simulink_demo.m` |
| `results/` | outputs of `run_paper_study` and `run_wire_study`, produced with GNU Octave 8.4 |

## Simulink

```matlab
setup_paths
build_cable_simulink_model                        % creates cable_thermal_sim.slx
out = sim('cable_thermal_sim');                   % T_probes: centre, cond-ins, ins-jacket, jacket-air
run_simulink_demo                                 % builds, simulates and compares with ct_simulate
```

The model runs as fixed-step discrete with a sample time equal to the model time step. It has
these blocks:

* **Clock, Bias of dt/2, 1-D Lookup Table (flat):** the Table 2 current steps, evaluated at
  mid-step like `ct_simulate`.
* **Constant:** the air temperature.
* **MATLAB Function block:** one FDM/TDMA step through `ct_simulink_step`, which calls the same
  `ct_coefficients` and `ct_tdma` functions.
* **Unit Delay** (initial condition T0): holds the node-temperature state.
* **Selector, Scope and To Workspace:** the outputs.

To simulate another cable or other options, pass them to the builder:

```matlab
build_cable_simulink_model(ct_cable_preset('cable150'), 'cable150_sim', ct_options('air_speed', 0.2))
```

**Important:** there is no MATLAB or Simulink in the environment where this code was written.
The MATLAB functions and the three scripts were run and tested with GNU Octave 8.4. For the
Simulink part, two checks were possible without Simulink:

* `run_tests` emulates the discrete loop of the Simulink model (Unit Delay, flat lookup of the
  current at mid-step, and the MATLAB Function block's `ct_simulink_step`). The result is
  identical to `ct_simulate` (difference 0.0 K over 860 steps).
* `build_cable_simulink_model.m` was dry-run with stub versions of the Simulink functions, which
  confirmed that every block, parameter, line and the MATLAB Function script are constructed
  correctly.

The real Simulink calls have not been executed. They use standard documented APIs
(`new_system`, `add_block`, `add_line`, `sfroot` / `Stateflow.EMChart.Script`), but treat the
builder and `run_simulink_demo.m` as untested until they have run once in your MATLAB.

## Verification (`run_tests`)

* TDMA vs a dense solve.
* Mesh areas, perimeters and conductor mass (Eqs. 6–24).
* Discrete energy balance (about 1e-9 W/m).
* Steady state vs the analytic cylinder solution.
* Adiabatic short circuit vs IEC 60949 (K = 226).
* Eq. 25 corrected = IEEE 738 q_cn, and the Eq. 26 units.
* Churchill–Chu against a hand calculation.
* The IEC 60287 iteration.
* **Parity with Python:** both test sequences agree within 5e-5 K, the CSV rounding; the
  ampacities of 10 database wires agree within 1e-5 relative.
* **Simulink step emulation** = `ct_simulate`.

There are 19 checks in total, all PASS in Octave 8.4 (`results/run_tests_octave.log`).

## Results (GNU Octave 8.4)

Same numbers as the Python version:

* **70 mm²:** end of the test 44.77 / 43.47 / 40.12 °C, against the paper's 39.31 / 37.52 / 35.21 °C.
* **150 mm²:** 38.68 / 37.25 / 34.50 °C, against the paper's 39.54 / 37.85 / 37.00 °C.
* **Ampacity at 90 °C:** 338.2 A and 548.5 A (IEC 60287: 354.3 A and 561.0 A).
* **Short circuit, 9.27 kA for 1 s:** peak 233.8 °C (IEC 60949: 250.0 °C adiabatic, 231.9 °C non-adiabatic).
* **Wire campaign:** see `results/wires/comparisons.txt` and `wire_results_matlab.csv`.
