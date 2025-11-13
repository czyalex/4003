#!/usr/bin/env python3
# 汇总最新 formal_* 下每个药的 Top-10 HPO，并映射到 MedDRA PT，便于放 PPT
import csv
from pathlib import Path

ROOT=Path(".")
MAP = ROOT/"data/interim/mappings/hpo_meddra_map.tsv"
OUTFOLDER = Path("reports")/("midterm_"+__import__("datetime").datetime.now().strftime("%Y-%m-%d_%H%M"))
OUT = OUTFOLDER/"top10_preview.csv"

def latest_formal():
    dirs = sorted((ROOT/"reports").glob("formal_*"), key=lambda p: p.stat().st_mtime, reverse=True)
    return dirs[0] if dirs else None

def main():
    formal = latest_formal()
    assert formal, "No formal_* dir found."
    h2pt = {}
    if MAP.exists():
        for r in csv.DictReader(open(MAP,encoding="utf-8-sig"), delimiter="\t"):
            h=(r.get("hpo_id") or "").strip().upper()
            pt=(r.get("meddra_pt") or "").strip()
            if h and pt: h2pt[h]=pt
    rows=[]
    for f in sorted(formal.glob("CHEMBL*_ranked_hpo.csv")):
        drug = f.stem.split("_")[0]
        with open(f,encoding="utf-8-sig") as fin:
            rdr = csv.DictReader(fin)
            for i,r in enumerate(rdr, start=1):
                if i>10: break
                h=r["hpo_id"].upper(); sc=r["score"]
                rows.append({"drug":drug,"rank":i,"hpo_id":h,"pt":h2pt.get(h,""),"score":sc})
    OUTFOLDER.mkdir(parents=True, exist_ok=True)
    with open(OUT,"w",encoding="utf-8",newline="") as fo:
        w=csv.DictWriter(fo, fieldnames=["drug","rank","hpo_id","pt","score"])
        w.writeheader(); w.writerows(rows)
    print("Wrote", OUT)

if __name__=="__main__":
    main()
