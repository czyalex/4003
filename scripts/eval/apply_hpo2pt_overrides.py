#!/usr/bin/env python3
import csv
from pathlib import Path

BASE = Path("data/interim/mappings/hpo_meddra_map.tsv")
OVR  = Path("data/interim/mappings/hpo_meddra_overrides.tsv")

def load_map(p):
    m={}
    if p.exists():
        for r in csv.DictReader(open(p,encoding="utf-8-sig"), delimiter="\t"):
            h=(r.get("hpo_id") or "").strip().upper()
            pt=(r.get("meddra_pt") or "").strip()
            if h and pt: m[h]=pt
    return m

def write_map(p, m):
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p,"w",encoding="utf-8",newline="") as f:
        w=csv.writer(f, delimiter="\t")
        w.writerow(["hpo_id","meddra_pt"])
        for h,pt in sorted(m.items()):
            w.writerow([h,pt])

def main():
    base=load_map(BASE)
    ovr =load_map(OVR)
    base.update(ovr)   # 覆盖
    write_map(BASE, base)
    print(f"Applied {len(ovr)} overrides. Total entries: {len(base)}")

if __name__=="__main__":
    main()
