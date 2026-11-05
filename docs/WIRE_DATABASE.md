# Wire database and test campaign

The cable model of Riba & Llauradó (Materials 2022, 15, 6814) is applied to a database of real
manufacturer wires from **36 AWG to 1/0 AWG** with nine insulation families, and the results are
compared with manufacturer current ratings. Results: [`results/wires/RESULTS_WIRES.md`](../results/wires/RESULTS_WIRES.md).

```bash
python -m cable_thermal wire-list --awg 22 --insulation silicone
```

```bash
python -m cable_thermal wire --part 3051 --current 8 --duration 300
```

```bash
python -m cable_thermal wire-study --out results/wires
```

## 1. Sources

| Data set | What | Where it comes from |
|---|---|---|
| `cable_thermal/data/wires/alphawire_hookup.json` | **337 Alpha Wire hook-up wire parts**: AWG, stranding, conductor plating and diameter, insulation material(s) and wall, OD, jacket, temperature range, voltage, DC resistance at 20 °C, UL AWM styles, current rating where published | Alpha Wire customer specification PDFs, `https://www.alphawire.com/disteAPI/SpecPDF/DownloadProductSpecPdf?productPartNumber=<part>` (part list from the site map), retrieved 2026-11-05. Rebuild with `tools/fetch_alphawire.py`. Families: Premium hook-up wire, EcoWire / EcoWire+ (mPPE), ThermoThin, PV wire. |
| `cable_thermal/data/wires/reference_ratings.json` | Wires **with published current ratings** | Habia Cable type E PTFE equipment wire (32–12 AWG, 25 types, ratings @ 40 °C + ambient correction table), [doc EQ-eE-04](https://www.triopak.fi/files/product/Habia/512003207.pdf); RS PRO PTFE type A (26/24/22 AWG, 40 K-rise ratings), [1](https://docs.rs-online.com/a078/0900766b8135dba2.pdf) [2](https://docs.rs-online.com/bff6/0900766b8135dba3.pdf) [3](https://docs.rs-online.com/68ed/0900766b8135dba4.pdf); General Cable Carol Super Vu-Tron 90 °C welding cable (6 AWG–4/0 AWG), [data sheet](https://www.datasheets.com/general-cable/01761.38.04/datasheet.pdf). Parsed with `tools/parse_reference_ratings.py`. |
| `cable_thermal/data/insulation_materials.json` | k, density, c_p of each insulation material, with the source of every number | see §3 |

Only facts (dimensions, ratings) are stored; the manufacturers' documents are not redistributed.

Coverage (Alpha Wire): PVC 123 parts (32–1/0 AWG, 80/90/105 °C), silicone rubber 73 (22–1 AWG,
150/200 °C, up to 40 kV), PTFE 58 (32–6 AWG, 200/250 °C), mPPE 24 (28–10 AWG), XLPE 20 (22–4 AWG,
incl. PV wire with PVC jacket), ETFE 11, ThermoThin fluoroplastic 11 (36–16 AWG, 300 °C), PVDF 4,
EPR 3, mica-glass + PTFE-fiberglass 10 (not simulated: no thermal data).

## 2. From a datasheet to a model

| Item | Choice |
|---|---|
| Conductor diameter D1 | datasheet value |
| Metal area | strands × π/4 · d², strand diameter from the AWG definition d = 0.127 mm · 92^((36−n)/39) (ASTM B258) or the stated inch value |
| r20 | datasheet DC resistance (Ω/1000 ft at 20 °C, nominal); 2 parts without DCR: ρ_Cu/area |
| α | 0.00393 1/K (annealed copper at 20 °C) |
| Conductor heat capacity | 8960 kg/m³ × metal area × c_p,Cu(T) (paper curve); plating neglected |
| Single insulation layer | two equal sub-layers of the same material (paper scheme with B = C, mathematically a single layer) |
| Insulation + jacket | layer B and layer C of the paper |
| Outer diameter | outermost wall adjusted so that the model OD equals the datasheet OD (flagged when the change is > 5 %) |
| Surface | horizontal wire, still air, ε = 0.85 (paper), natural convection **Churchill & Chu** (see below), sea level |
| Mesh / time step | 9 nodes per material (paper); transient Δt scaled to the wire's time constant |

**Convection law.** The paper's Eq. 25 (IEEE 738, h ∝ D^−0.25 ΔT^0.25) is a power-law fit for
overhead conductors (Ra ≈ 10²–10⁴⁺). Hook-up wires reach Ra < 1, where it under-predicts h; the
Churchill & Chu (1975) correlation is valid for 10⁻⁵ ≤ Ra ≤ 10¹². Ratio of the ampacity with Eq. 25
to that with Churchill–Chu over the database: median 0.84 for OD < 1 mm, 0.91 for 1–2 mm, 0.97 for
2–5 mm, 1.01 above 5 mm. The RS PRO ratings (defined by a 40 K rise) are reproduced with
Churchill–Chu (35–43 K) but not with Eq. 25 (52–58 K), so Churchill–Chu is the default for wires
(`--convection` switches it).

## 3. Insulation materials

| Key | k [W/(m K)] | ρ [kg/m³] | c_p [J/(kg K)] | Sources |
|---|---|---|---|---|
| PVC | 0.142 → 0.113 (20 → 75 °C) | 1400* | 950 → 1180 | paper Figs. 2c/3c; *density assumed |
| XLPE | 0.31–0.33 (20–75 °C) | 920* | 1760 → 3150 | paper Figs. 2b/3b; *density assumed |
| PTFE | 0.25 | 2170 | 1000 | Fluorotherm material comparison; fluororesin characteristics table (JIS K 7123) |
| FEP / PFA | 0.20 / 0.19 | 2150 | 1200 / 1000 | same |
| ETFE | 0.24 | 1740 | 2000 | same |
| PVDF | 0.18 | 1780 | 1200 | Arkema Kynar 740 TDS (k 0.17–0.19, ρ 1.77–1.79); fluororesin table (c_p) |
| ThermoThin "fluoroplastic" | 0.22 | 2150 | 1100 | type not disclosed by Alpha – mean of PTFE/FEP/PFA/ETFE |
| Silicone rubber | 0.20 | 1200 | 1256 | Shin-Etsu silicone rubber catalog (≈0.2); TYM silicone sheet (SG 1.2); Momentive RTV615 (0.3 cal/g °C) |
| mPPE | 0.27 | 1110 | 1200 | SABIC NORYL V0150B (ISO 8302); Ensinger TECANYL GF30 (c_p, ISO 22007-4) |
| EPR | 0.286 | 1300 | 1538 | IEC 60287-2-1 (ρ_T = 3.5 K m/W), IEC 60853-2 (ρc = 2.0 MJ/(m³ K)) |

Mappings: irradiated/semi-rigid PVC → PVC; PTFE tape → PTFE; mPPE-TPE → mPPE; the 0.005" nylon
over-layer of THHN → PVC properties (no nylon data sourced); mica-glass tape and PTFE-impregnated
fiberglass → not simulated. Above 75 °C the paper's PVC/XLPE/Cu curves are held constant. The
steady-state ampacity depends only on k; k ±30 % changes it by −6 … +4 % (larger wires most).

## 4. What the test campaign does

For every simulatable part (327): ampacity with the conductor at the wire's rated temperature in
30 °C still air (Churchill–Chu and Eq. 25), ampacity at 105 °C, current for a 40 K rise, thermal
time constant (63 %), insulation temperature drop → `results/wires/wire_results.csv`.

Comparisons with published ratings (`results/wires/RESULTS_WIRES.md` §3, `fig23`, `fig24`):

| Data set | Criterion of the rating | Model result |
|---|---|---|
| RS PRO PTFE type A, 26–22 AWG | 40 K rise, single wire in free air | rise at the rating 34.7 / 39.5 / 42.9 K |
| Alpha Wire PV wire 8/6/4 AWG (XLPE + PVC jacket) | 90 °C wire, current @ 30 °C | conductor 77.3–77.6 °C at the rating; rating = 0.89 × model 90 °C ampacity |
| Super Vu-Tron welding cable 6 AWG–4/0 AWG | 90 °C, NEC 400.5(B), 30 °C | conductor 75.6–78.0 °C at the rating; rating = 0.87–0.89 × model ampacity (conductor diameter estimated, packing 0.75–0.85 changes it by ±1 %) |
| Habia PTFE type E, 32–12 AWG | "Amps @ 40 °C", 200 °C wire, basis not stated | conductor 112–191 °C at the rating (median 148 °C); rating = 0.69–0.97 × model 200 °C ampacity. Habia's ambient correction is steeper than the model's for a 200 °C conductor (0.28 vs 0.36 at 180 °C), consistent with a conductor limit below 200 °C |

## 5. Limitations

* Single horizontal wire in still air, no bundling, no enclosure, no solar load.
* Emissivity 0.85 for all colours and materials (±0.1 → about ±2 % ampacity).
* Datasheet values are nominal; tolerance on wall and DCR is not propagated.
* Two EcoWire+ parts (6830, 6831) list a DC resistance that implies a resistivity below that of
  copper for the stated stranding – the datasheet value is used and flagged.
* Material properties are typical resin values, not those of the manufacturer's compound.
