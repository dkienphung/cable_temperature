"""Rebuild the Alpha Wire hook-up wire database from the manufacturer's web site.

Steps (network access to www.alphawire.com and `pdftotext` from poppler-utils required):
  1. read https://www.alphawire.com/sitemap.xml and collect the wire part numbers per family;
  2. download each customer specification PDF from
     https://www.alphawire.com/disteAPI/SpecPDF/DownloadProductSpecPdf?productPartNumber=<part>
     (one request every 0.4 s);
  3. convert to text with `pdftotext -layout`;
  4. parse with tools/parse_alphawire.py and write cable_thermal/data/wires/alphawire_hookup.json
     (facts only - the PDFs are not stored in the repository).

Usage:  python tools/fetch_alphawire.py <work_dir>
"""
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120 Safari/537.36"
FAMILIES = ["ecogen/ecowire", "ecogen/ecowireplus", "hook-up-wire/thermothin", "hook-up-wire/premium", "hook-up-wire/pv"]


def curl(url, out):
    r = subprocess.run(["curl", "-sSL", "-A", UA, "--max-time", "60", "-o", out, "-w", "%{http_code} %{content_type}", url],
                       capture_output=True, text=True)
    return r.stdout.strip()


def main(work):
    os.makedirs(os.path.join(work, "pdf"), exist_ok=True)
    os.makedirs(os.path.join(work, "txt"), exist_ok=True)
    sm = os.path.join(work, "sitemap.xml")
    curl("https://www.alphawire.com/sitemap.xml", sm)
    locs = re.findall(r"<loc>([^<]+)</loc>", open(sm, encoding="utf-8").read())
    fam = {}
    for u in locs:
        p = u.split("/")
        if len(p) >= 8 and p[3].lower() == "products" and p[4].lower() == "wire":
            g = "/".join(p[5:-1]).lower()
            if g in FAMILIES:
                fam[p[-1].upper().replace("_", "/")] = [g]
    for part in sorted(fam):
        fn = os.path.join(work, "pdf", part.replace("/", "_") + ".pdf")
        if not (os.path.exists(fn) and os.path.getsize(fn) > 1000):
            st = curl("https://www.alphawire.com/disteAPI/SpecPDF/DownloadProductSpecPdf?productPartNumber=" + part, fn)
            fam[part].append(st)
            time.sleep(0.4)
        subprocess.run(["pdftotext", "-layout", fn, os.path.join(work, "txt", part.replace("/", "_") + ".txt")])
    json.dump(fam, open(os.path.join(work, "family_map.json"), "w"), indent=1)
    raw = os.path.join(work, "aw_raw.json")
    subprocess.run([sys.executable, os.path.join(HERE, "parse_alphawire.py"), os.path.join(work, "txt"),
                    os.path.join(work, "family_map.json"), raw], check=True)
    recs = [{k: v for k, v in r.items() if not k.startswith("raw_")} for r in json.load(open(raw))]
    date = time.strftime("%Y-%m-%d")
    for r in recs:
        r["retrieved"] = date
    out = os.path.join(HERE, "..", "cable_thermal", "data", "wires", "alphawire_hookup.json")
    json.dump({"_about": "Alpha Wire hook-up wire customer specifications, parsed by tools/parse_alphawire.py",
               "manufacturer": "Alpha Wire (Richmond, IN, USA)", "retrieved": date, "parts": recs},
              open(out, "w"), indent=1, ensure_ascii=False)
    print(len(recs), "parts written to", out)


if __name__ == "__main__":
    main(sys.argv[1])
