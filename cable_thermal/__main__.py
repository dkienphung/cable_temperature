"""Command-line application.

Examples
--------
python -m cable_thermal reproduce --out results
python -m cable_thermal simulate --cable cable70                      # paper test sequence (Table 2)
python -m cable_thermal simulate --cable cable70 --steps 150:3600,250:1800 --ambient 25 --plot run.png --csv run.csv
python -m cable_thermal simulate --config examples/overload.json
python -m cable_thermal steady --cable cable150 --current 400 --ambient 30
python -m cable_thermal ampacity --cable cable70 --tmax 90 --ambient 30
python -m cable_thermal simulate --cable cable70 --air-speed 0.2              # indoor draught 0.2 m/s
python -m cable_thermal ampacity --cable cable70 --elevation 0 --air-speed 0.5 --air-angle 45
python -m cable_thermal wire-list --awg 22                                  # wire database (Alpha Wire)
python -m cable_thermal wire --part 3051 --ambient 30                       # ampacity + temperature of a wire
python -m cable_thermal wire --part 391250 --current 8 --duration 600 --plot w.png
python -m cable_thermal wire-study --out results/wires                     # whole-database test campaign
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from .cables import get_cable
from .model import CableThermalModel, ModelOptions


def _model(args, cable=None):
    cable = cable or get_cable(args.cable)
    opts = ModelOptions(nodes_per_material=tuple([args.nodes] * 3), dt=args.dt, convection_model=args.convection,
                        elevation=args.elevation, air_speed=args.air_speed, air_flow_angle=args.air_angle)
    return CableThermalModel(cable, opts)


def _parse_steps(text):
    steps = []
    for part in text.split(","):
        I, d = part.split(":")
        steps.append((float(I), float(d)))
    return steps


def cmd_simulate(args):
    cable = get_cable(args.cable)
    current, ambient = cable.test_steps, args.ambient
    t_end = None
    if args.config:
        cfg = json.load(open(args.config))
        cable = get_cable(cfg.get("cable", args.cable))
        if "overrides" in cfg:
            cable = cable.with_(**cfg["overrides"])
        current = [tuple(s) for s in cfg["current_steps"]]
        ambient = cfg.get("ambient", ambient)
        if isinstance(ambient, list):
            ambient = [tuple(s) for s in ambient]
        t_end = cfg.get("t_end")
        for k, v in cfg.get("model", {}).items():
            setattr(args, k, v)
    elif args.steps:
        current = _parse_steps(args.steps)
    m = _model(args, cable)
    r = m.simulate(current, ambient=ambient, t_end=t_end)
    print(f"cable: {cable.designation} | H = {m.elevation:g} m, air speed = {m.opt.air_speed:g} m/s")
    print(f"{'t [s]':>8} {'I [A]':>7} {'T_air':>6} {'T_centre':>9} {'T_c-ins':>8} {'T_ins-jac':>9} {'T_jac-air':>9}")
    step = max(1, len(r.t) // 20)
    for j in sorted(set(range(0, len(r.t), step)) | {len(r.t) - 1}):
        print(f"{r.t[j]:8.0f} {r.current[j]:7.1f} {r.ambient[j]:6.1f} {r.T_center[j]:9.2f} {r.T_conductor_insulation[j]:8.2f} "
              f"{r.T_insulation_jacket[j]:9.2f} {r.T_jacket_air[j]:9.2f}")
    if args.csv:
        r.to_csv(args.csv)
        print("csv written to", args.csv)
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for name, lab, col in [("conductor_insulation", "Conductor-insulation", "#eb6834"),
                               ("insulation_jacket", "Insulation-jacket", "#2a78d6"),
                               ("jacket_air", "Jacket-air", "#6b6a65")]:
            ax.plot(r.t, r.probe(name), color=col, label=lab)
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Temperature [°C]")
        ax.grid(alpha=0.3)
        ax.legend()
        ax.set_title(cable.designation, loc="left", fontsize=10)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=200)
        print("plot written to", args.plot)


def cmd_steady(args):
    m = _model(args)
    res = m.steady_heat_terms(args.current, args.ambient)
    T = res["T"]
    print(f"I = {args.current:g} A, T_air = {args.ambient if args.ambient is not None else m.cable.ambient_temperature:g} °C, "
          f"H = {m.elevation:g} m, air speed = {m.opt.air_speed:g} m/s")
    for k, i in m.probes.items():
        print(f"  {k:22s} {T[i]:8.2f} °C")
    print(f"  Joule {res['joule']:.3f} W/m = convection {res['p_conv']:.3f} + radiation {res['p_rad']:.3f} W/m "
          f"(h_conv = {res['h_conv']:.2f} W/m²K)")


def cmd_ampacity(args):
    m = _model(args)
    I = m.ampacity(args.tmax, args.ambient)
    print(f"steady-state current giving {args.tmax:g} °C at the conductor: {I:.1f} A")


def cmd_wire_list(args):
    from .wires import all_alpha_wires, awg_sort_key

    rows = [w for w in all_alpha_wires()
            if (args.awg is None or w.awg == args.awg)
            and (args.insulation is None or args.insulation.lower() in w.insulation_label.lower())]
    rows.sort(key=lambda w: (-awg_sort_key(w.awg), w.insulation_label, w.part_number))
    print(f"{'part':>10} {'AWG':>4} {'insulation':32s} {'Tmax':>5} {'D [mm]':>7} {'OD [mm]':>7} {'R20 [ohm/km]':>12}")
    for w in rows:
        print(f"{w.part_number:>10} {w.awg:>4} {w.insulation_label[:32]:32s} {str(w.temp_max_C):>5} "
              f"{w.conductor_diameter * 1e3:7.3f} {w.outer_diameter * 1e3:7.3f} {w.r20 * 1e3:12.2f}")
    print(f"{len(rows)} parts")


def cmd_wire(args):
    from .wires import find_part

    w = find_part(args.part)
    print(f"{w.source} {w.part_number}: {w.awg} AWG, {w.insulation_label}, rated {w.temp_max_C} °C, "
          f"D = {w.conductor_diameter * 1e3:.3f} mm, OD = {w.outer_diameter * 1e3:.3f} mm, "
          f"R20 = {w.r20 * 1e3:.2f} ohm/km ({w.r20_source})")
    for n in w.notes:
        print("  note:", n)
    if not w.simulatable:
        print("  not simulated: no thermal data for this insulation")
        return
    c = w.to_cable(args.ambient, elevation=args.elevation or 0.0)
    opts = ModelOptions(nodes_per_material=tuple([args.nodes] * 3), dt=args.dt, convection_model=args.convection,
                        elevation=args.elevation, air_speed=args.air_speed, air_flow_angle=args.air_angle)
    m = CableThermalModel(c, opts)
    tmax = args.tmax if args.tmax is not None else float(w.temp_max_C)
    print(f"  ampacity (conductor {tmax:g} °C, air {args.ambient:g} °C, {args.convection}): {m.ampacity(tmax, args.ambient):.2f} A")
    if args.current is not None:
        T = m.steady_state(args.current, args.ambient)
        print(f"  steady state at {args.current:g} A: conductor {T[0]:.2f} °C, surface {T[-1]:.2f} °C")
        if args.duration:
            r = m.simulate(args.current, ambient=args.ambient, t_end=args.duration, dt=min(args.dt, args.duration / 600))
            for tt in np.linspace(0, args.duration, 7):
                print(f"    t = {tt:8.1f} s  conductor {np.interp(tt, r.t, r.T_center):7.2f} °C  surface {np.interp(tt, r.t, r.T_jacket_air):7.2f} °C")
            if args.plot:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(7, 4))
                ax.plot(r.t, r.T_center, color="#eb6834", label="conductor")
                ax.plot(r.t, r.T_jacket_air, color="#6b6a65", label="surface")
                ax.set_xlabel("Time [s]")
                ax.set_ylabel("Temperature [°C]")
                ax.set_title(f"{w.part_number} {w.awg} AWG {w.insulation_label}, {args.current:g} A", loc="left", fontsize=10)
                ax.legend()
                ax.grid(alpha=0.3)
                fig.tight_layout()
                fig.savefig(args.plot, dpi=200)
                print("  plot written to", args.plot)


def cmd_wire_study(args):
    from .wire_study import run_wire_study

    rows, _ = run_wire_study(args.out)
    print(f"{len(rows)} parts processed; report {args.out}/RESULTS_WIRES.md, table {args.out}/wire_results.csv")


def cmd_reproduce(args):
    from .reproduce import run_all

    s = run_all(args.out)
    print(json.dumps({k: s[k] for k in ("comparison_with_paper", "diagnostic_calibration", "steady_state")}, indent=1, default=float)[:4000])
    print(f"figures in {args.out}/figures, report {args.out}/RESULTS.md")


def main(argv=None):
    p = argparse.ArgumentParser(prog="cable_thermal", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--cable", default="cable70", help="cable70 | cable150")
        sp.add_argument("--nodes", type=int, default=9, help="nodes per material (paper: 9)")
        sp.add_argument("--dt", type=float, default=10.0, help="time step [s] (paper: 10)")
        sp.add_argument("--convection", default="ieee738", choices=["ieee738", "paper", "churchill_chu"],
                        help="ieee738 = corrected Eq. 25 (rho^0.5); paper = Eq. 25 as printed (rho^0.25); "
                             "churchill_chu = horizontal-cylinder correlation valid for thin wires")
        sp.add_argument("--elevation", type=float, default=None,
                        help="site elevation H in Eq. 26 [m] (default: cable preset, 280 m for the paper's lab)")
        sp.add_argument("--air-speed", dest="air_speed", type=float, default=0.0,
                        help="indoor air speed [m/s]; 0 = still air as in the paper, >0 adds IEEE 738 forced convection")
        sp.add_argument("--air-angle", dest="air_angle", type=float, default=90.0,
                        help="angle between air flow and cable axis [deg] (default 90 = perpendicular)")

    s = sub.add_parser("simulate", help="transient / dynamic simulation")
    common(s)
    s.add_argument("--steps", help="current steps 'I:duration,I:duration' (default: paper test sequence)")
    s.add_argument("--ambient", type=float, default=None)
    s.add_argument("--config", help="JSON scenario file")
    s.add_argument("--csv")
    s.add_argument("--plot")
    s.set_defaults(func=cmd_simulate)

    s = sub.add_parser("steady", help="steady-state temperatures for a current")
    common(s)
    s.add_argument("--current", type=float, required=True)
    s.add_argument("--ambient", type=float, default=None)
    s.set_defaults(func=cmd_steady)

    s = sub.add_parser("ampacity", help="current for a maximum conductor temperature")
    common(s)
    s.add_argument("--tmax", type=float, default=90.0)
    s.add_argument("--ambient", type=float, default=None)
    s.set_defaults(func=cmd_ampacity)

    s = sub.add_parser("wire-list", help="list the wire database")
    s.add_argument("--awg")
    s.add_argument("--insulation")
    s.set_defaults(func=cmd_wire_list)

    s = sub.add_parser("wire", help="ampacity / temperature of a database wire")
    common(s)
    s.set_defaults(convection="churchill_chu")
    s.add_argument("--part", required=True)
    s.add_argument("--ambient", type=float, default=30.0)
    s.add_argument("--tmax", type=float, default=None, help="conductor limit (default: wire rating)")
    s.add_argument("--current", type=float, default=None)
    s.add_argument("--duration", type=float, default=None)
    s.add_argument("--plot")
    s.set_defaults(func=cmd_wire)

    s = sub.add_parser("wire-study", help="run the model on the whole wire database and compare with ratings")
    s.add_argument("--out", default="results/wires")
    s.set_defaults(func=cmd_wire_study)

    s = sub.add_parser("reproduce", help="reproduce all paper figures + verification report")
    s.add_argument("--out", default="results")
    s.set_defaults(func=cmd_reproduce)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
