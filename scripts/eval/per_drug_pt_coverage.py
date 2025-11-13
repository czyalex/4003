#!/usr/bin/env python3
# Report per-drug coverage: among top-K HPO, how many map to a MedDRA PT.
import csv
from pathlib import Path
from datetime import date

K = 10
MAP = Path("data/interim/mappings/hpo_meddra_map.tsv")
OUTDIR = Path("reports")/f"formal_{date.today()}"; OUTDIR.mkdir(parents=True, exist_ok=True)

def read_map():
    d = {}
    with MAP.open("r",encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            d[row["hpo_id"]] = row["meddra_pt"]
    return d

def main():
    h2pt = read_map()
    rows = []
    for p in OUTDIR.glob("CHEMBL*_ranked_hpo.csv"):
        drug = p.stem.split("_")[0]
        top = []
        with p.open("r",encoding="utf-8") as f:
            rd = csv.DictReader(f)
            for i,row in enumerate(rd, start=1):
                top.append(row["hpo_id"])
                if i>=K: break
        cov = sum(1 for h in top if h in h2pt)
        rows.append({"drug":drug,"K":K,"mapped":cov,"coverage":f"{cov/K:.3f}"})

    with (OUTDIR/"pt_coverage.csv").open("w",encoding="utf-8",newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=["drug","K","mapped","coverage"])
        w.writeheader(); w.writerows(rows)
    print(f"Wrote {(OUTDIR/'pt_coverage.csv')}")

if __name__=="__main__":
    main()
