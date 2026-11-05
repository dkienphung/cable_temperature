"""Parse Alpha Wire customer-specification sheets (pdftotext -layout output) into records.

Usage:  python tools/parse_alphawire.py <txt_dir> <family_map.json> <out.json>

Each record keeps the raw datasheet lines it was parsed from, so every number can be traced
back to the manufacturer's document.
"""
import json
import os
import re
import sys

IN = 25.4e-3  # m per inch


def _num(s):
    return float(s) if s is not None else None


def parse(text, part):
    rec = {"part_number": part, "manufacturer": "Alpha Wire"}
    lines = text.splitlines()
    # ---------------------------------------------------------------- conductor
    cond = next((l for l in lines if re.match(r"\s*a\)\s*Conductor", l)), None)
    if cond is None:
        rec["parse_error"] = "no conductor line"
        return rec
    rec["raw_conductor"] = re.sub(r"\s+", " ", cond.strip())
    m = re.search(r"Conductor\s+(\d+(?:/0)?)\s*\((.*)\)\s*AWG\s+(.*?)\s{2,}([\d.]+)\s*$", cond)
    if not m:
        m = re.search(r"Conductor\s+(\d+(?:/0)?)\s*\((.*)\)\s*AWG\s+(.*?)\s+([\d.]+)\s*$", cond)
    if not m:
        rec["parse_error"] = "conductor format"
        return rec
    awg, strand, metal, d = m.groups()
    rec["awg"] = awg
    rec["stranding"] = strand.strip()
    rec["conductor_material"] = re.sub(r"\s+", " ", metal.strip()).replace("Platted", "Plated")
    rec["conductor_diameter_in"] = float(d)
    m2 = re.match(r"^\s*(\d+)\s*/\s*(\d+)\s*$", strand)
    if strand.strip().upper() == "SOLID":
        rec["strands"], rec["strand_awg"] = 1, int(awg) if awg.isdigit() else None
    elif m2:
        rec["strands"], rec["strand_awg"] = int(m2.group(1)), int(m2.group(2))
    # ---------------------------------------------------------------- insulation (+ continuation)
    i_ins = next((k for k, l in enumerate(lines) if re.match(r"\s*b\)\s*Insulation", l)), None)
    if i_ins is not None:
        ins = lines[i_ins]
        rec["raw_insulation"] = re.sub(r"\s+", " ", ins.strip())
        mi = re.search(r'Insulation\s+([\d.]+)"\s*Wall,\s*Nom\.?\s*,?\s*(.*?)\s{2,}([\d.]+)\s*(?:\+/-\s*([\d.]+))?', ins)
        if mi:
            wall, mat, od, tol = mi.groups()
            layers = [{"material": mat.strip(), "wall_in": float(wall)}]
            # e.g. "PVC/ 0.005" Wall NYLON" -> second layer
            m3 = re.match(r'(.*?)/\s*([\d.]+)"\s*Wall\s*,?\s*(?:Nom\.?)?\s*(.*)', mat)
            if m3:
                layers = [{"material": m3.group(1).strip(), "wall_in": float(wall)},
                          {"material": m3.group(3).strip(), "wall_in": float(m3.group(2))}]
            rec["insulation_layers"] = layers
            rec["insulation_od_in"] = float(od)
            rec["insulation_od_tol_in"] = _num(tol)
    # ---------------------------------------------------------------- jacket (component 2)
    for k, l in enumerate(lines):
        mj = re.search(r'^\s*2\)\s*Jacket\s+([\d.]+)"\s*Wall,\s*Nom\.?\s*,?\s*(.*?)\s{2,}([\d.]+)', l)
        if mj:
            rec["jacket"] = {"material": mj.group(2).strip(), "wall_in": float(mj.group(1)), "od_in": float(mj.group(3))}
            rec["raw_jacket"] = re.sub(r"\s+", " ", l.strip())
    if any("Separator" in l for l in lines):
        rec["separator"] = next(re.sub(r"\s+", " ", l.split("Separator", 1)[1].strip()) for l in lines if "Separator" in l)
    # ---------------------------------------------------------------- properties
    t = re.search(r"Temperature Range\s+(-?\d+)\s*to\s*(-?\d+)\s*°C", text)
    if t:
        rec["temp_min_C"], rec["temp_max_C"] = int(t.group(1)), int(t.group(2))
    v = re.search(r"Voltage Rating\s+([\d,]+)\s*V", text)
    if v:
        rec["voltage_V"] = int(v.group(1).replace(",", ""))
    r = re.search(r"Conductor DCR\s+([\d.]+)\s*Ω/1000\s*ft\s*@\s*20°C", text)
    if r:
        rec["dcr_ohm_per_kft_20C"] = float(r.group(1))
    c = re.search(r"Current\s+(\d+)\s*amps per conductor\s*@\s*(\d+)°C", text)
    if c:
        rec["current_rating_A"], rec["current_rating_ambient_C"] = int(c.group(1)), int(c.group(2))
    styles = sorted(set(re.findall(r"(?:STYLES?|AWM)\s+(\d{4,5})", text)))
    if styles:
        rec["ul_awm_styles"] = styles
    return rec


def main(txt_dir, family_map, out):
    fam = json.load(open(family_map))
    recs = []
    for fn in sorted(os.listdir(txt_dir)):
        if not fn.endswith(".txt"):
            continue
        part = fn[:-4].replace("_", "/")
        rec = parse(open(os.path.join(txt_dir, fn), encoding="utf-8", errors="replace").read(), part)
        rec["family"] = fam.get(part, [None])[0]
        rec["source_url"] = ("https://www.alphawire.com/disteAPI/SpecPDF/DownloadProductSpecPdf?productPartNumber="
                             + part)
        recs.append(rec)
    json.dump(recs, open(out, "w"), indent=1, ensure_ascii=False)
    print(len(recs), "records,", sum("parse_error" in r for r in recs), "parse errors")


if __name__ == "__main__":
    main(*sys.argv[1:4])
