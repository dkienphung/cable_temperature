"""Parse third-party manufacturer datasheets that publish current ratings (pdftotext -layout text).

Usage: python tools/parse_reference_ratings.py <habia_E.txt> <rs_a1.txt> <rs_a2.txt> <rs_a3.txt> <supervutron.txt> <out.json>
"""
import json
import re
import sys

IN = 25.4e-3


def habia_e(text):
    rows = []
    for line in text.splitlines():
        m = re.match(r"\s*E\s+(\d{4})\s+(\d+)\s+([\d,]+)\s+(\d+)\s*x\s*([\d,]+)\s+([\d,]+)\s+([\d,]+)\s+([\d,]+)\s+±\s*([\d,]+)\s+([\d,]+)\s+(\d+)\s+(\S+)", line)
        if m:
            f = lambda s: float(s.replace(",", "."))
            rows.append({"type": "E " + m.group(1), "awg": m.group(2), "area_mm2": f(m.group(3)),
                         "strands": int(m.group(4)), "strand_diameter_mm": f(m.group(5)),
                         "conductor_diameter_mm": f(m.group(6)), "dcr_ohm_per_km_20C": f(m.group(7)),
                         "od_mm": f(m.group(8)), "od_tol_mm": f(m.group(9)), "mass_kg_per_km": f(m.group(10)),
                         "current_rating_A": int(m.group(11)), "order_ref": m.group(12)})
    m = re.search(r"°C\s+((?:\d+\s+){8}\d+)\s*\n.*?\n?\s*Factor\s+((?:[\d.]+\s+){8}[\d.]+)", text)
    derate = None
    if m:
        derate = dict(zip([int(x) for x in m.group(1).split()], [float(x) for x in m.group(2).split()]))
    return rows, derate


def rs_ptfe(text):
    m = re.search(r"\n\s*(\d{2})\s+([\d.]+)\s+([\d.]*)\s*(\d+/[\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", text)
    if not m:
        return None
    g = m.groups()
    od = g[2] if g[2] else None
    m_od = re.search(r"O/D Stranding.*?\n\s*([\d.]+)\s*\n", text, flags=re.S)
    cond_od = float(od) if od else (float(m_od.group(1)) if m_od else None)
    title = re.search(r"Equipment Wire \(Type A\)\s*\n\s*([\d/.]+mm SPC)", text)
    return {"product": "RS PRO PTFE Insulated Equipment Wire (Type A) " + (title.group(1) if title else ""),
            "awg": g[0], "area_mm2": float(g[1]), "conductor_diameter_mm": cond_od, "stranding": g[3],
            "insulation_thickness_mm_min": float(g[4]), "insulation_thickness_mm_max": float(g[5]),
            "od_mm_min": float(g[6]), "od_mm_max": float(g[7]), "dcr_ohm_per_km_20C_max": float(g[8]),
            "mass_kg_per_km": float(g[9]), "current_rating_A": float(g[10])}


def supervutron(text):
    rows = []
    for m in re.finditer(r"(0176\d)†?\s+(\d/0|\d)\s+(\d+)/34\s+([\d.]+)\s+([\d.]+)\s+(\d+)", text):
        rows.append({"catalog": m.group(1), "awg": m.group(2), "strands": int(m.group(3)), "strand_awg": 34,
                     "od_in": float(m.group(4)), "od_mm": float(m.group(5)), "net_weight_lb_per_kft": int(m.group(6))})
    amp = {}
    i = text.find("For 600 Volt In-Line Applications")
    seg = text[i:i + 1500]
    for a, b in re.findall(r"(\d/0|\b\d)\s+(\d{3})\b", seg):
        amp.setdefault(a, int(b))
    for r in rows:
        r["current_rating_A"] = amp.get(r["awg"])
    return rows


def main(habia, rs1, rs2, rs3, svt, out):
    hrows, derate = habia_e(open(habia).read())
    rs = [rs_ptfe(open(p).read()) for p in (rs1, rs2, rs3)]
    sv = supervutron(open(svt).read())
    data = {
        "habia_E_ptfe": {
            "manufacturer": "Habia Cable", "product": "Equipment Wire E, 600 V, +200 °C, NEMA HP-3, extruded PTFE, silver plated copper",
            "source_url": "https://www.triopak.fi/files/product/Habia/512003207.pdf", "document_ref": "EQ-eE-04 (2007-05-11)",
            "insulation": "PTFE", "conductor": "Silver Plated Copper", "temp_max_C": 200, "voltage_V": 600,
            "rating_basis": "'Amps (@ 40 °C)' = ambient air 40 °C; single wire assumed (multicore table gives 1.0 for 2 cores); conductor-temperature basis not stated",
            "rating_ambient_C": 40.0, "ambient_correction": derate, "rows": hrows},
        "rs_ptfe_type_a": {
            "manufacturer": "RS PRO (RS Components)", "insulation": "PTFE", "conductor": "Silver plated copper",
            "temp_max_C": 200, "voltage_V": 300, "standard": "BS 3G 210:1996",
            "rating_basis": "temperature rise of 40 °C, single wire in free air (ambient not stated)",
            "rating_temperature_rise_K": 40.0,
            "source_urls": ["https://docs.rs-online.com/a078/0900766b8135dba2.pdf", "https://docs.rs-online.com/bff6/0900766b8135dba3.pdf",
                            "https://docs.rs-online.com/68ed/0900766b8135dba4.pdf"],
            "rows": [r for r in rs if r]},
        "general_cable_super_vutron": {
            "manufacturer": "General Cable / Prysmian", "product": "Carol Super Vu-Tron Welding Cable, 90 °C 600 V UL/CSA RHH/RHW, Class M (34 AWG strands)",
            "source_url": "https://www.datasheets.com/general-cable/01761.38.04/datasheet.pdf",
            "insulation": "Super Vu-Tron (proprietary rubber jacket; single layer over the conductor)", "temp_max_C": 90,
            "rating_basis": "'Suggested ampacities for 600 V in-line applications', single conductor, per ICEA S-19-81 / NEMA WC-3 Appendix J and NEC Table 400.5(B); NEC ampacity tables use 30 °C ambient",
            "rating_ambient_C": 30.0, "rows": sv},
    }
    json.dump(data, open(out, "w"), indent=1, ensure_ascii=False)
    print(len(hrows), "Habia rows; derate", derate, "| RS", [(r["awg"], r["current_rating_A"]) for r in data["rs_ptfe_type_a"]["rows"]],
          "| SVT", [(r["awg"], r["current_rating_A"], r["od_in"]) for r in sv])


if __name__ == "__main__":
    main(*sys.argv[1:7])
