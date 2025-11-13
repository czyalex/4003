#!/usr/bin/env python3
# Quick lineage report: counts, unique entities, and duplicates in key mapping tables.

import csv
from pathlib import Path
from datetime import date

BASE = Path("data/interim/mappings")
OUTDIR = Path("reports")/f"formal_{date.today()}"; OUTDIR.mkdir(parents=True, exist_ok=True)
OUT = OUTDIR/"lineage_summary.txt"

def load(p, delim=","):
    return list(csv.DictReader(p.open("r",encoding="utf-8"), delimiter=delim))

def dedup_count(rows, keys):
    seen=set(); dup=0
    for r in rows:
        k=tuple(r[k] for k in keys); dup += 1 if k in seen else 0; seen.add(k)
    return dup

def main():
    dt = load(BASE/"drug_targets.csv")
    pg = load(BASE/"protein_gene.csv")
    gh = load(BASE/"gene_hpo.csv")
    mp = load(BASE/"hpo_meddra_map.tsv", "\t")

    lines = []
    lines.append(f"Drug→Protein rows: {len(dt)} | dup pairs: {dedup_count(dt,['drug_id','protein_id'])}")
    lines.append(f"Protein→Gene rows: {len(pg)} | dup pairs: {dedup_count(pg,['protein_id','gene_id'])}")
    lines.append(f"Gene→HPO rows   : {len(gh)} | dup pairs: {dedup_count(gh,['gene_id','hpo_id'])}")
    lines.append(f"HPO→PT rows     : {len(mp)} | dup pairs: {dedup_count(mp,['hpo_id','meddra_pt'])}")

    OUT.write_text("\n".join(lines)+"\n", encoding="utf-8")
    print(f"Wrote {OUT}")

if __name__=="__main__":
    main()
