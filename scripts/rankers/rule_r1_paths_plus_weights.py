#!/usr/bin/env python3
# Sum over all Drug→Protein→Gene→HPO paths with action weights, then add β·CPIC and λ·PT prior.

import csv, math
from pathlib import Path
from collections import defaultdict
from datetime import date

BETA   = 0.6   # CPIC prior weight
LAMBDA = 0.3   # PT frequency prior weight

BASE = Path("data/interim/mappings")
DT = BASE/"drug_targets.csv"
PG = BASE/"protein_gene.csv"
GH = BASE/"gene_hpo.csv"
ACTIONS = BASE/"drug_target_actions.tsv"
PGX = BASE/"drug_gene_pgx.csv"
HPO2PT = BASE/"hpo_meddra_map.tsv"
PTPRIOR = BASE/"pt_prior.csv"

OUTDIR = Path("reports")/f"formal_{date.today()}"; OUTDIR.mkdir(parents=True, exist_ok=True)

def read_csv(p): return list(csv.DictReader(p.open("r",encoding="utf-8")))
def read_tsv(p): return list(csv.DictReader(p.open("r",encoding="utf-8"), delimiter="\t"))

def action_weight(a: str) -> float:
    if not a: return 1.0
    a = a.lower()
    # simple mapping; tune later
    if "inhib" in a or "block" in a or "antagon" in a: return 1.2
    if "agon" in a: return 1.1
    if "modulat" in a: return 1.05
    return 1.0

def main():
    dt = read_csv(DT)
    pg = read_csv(PG)
    gh = read_csv(GH)

    # maps
    drug2prot = defaultdict(list)
    for r in dt:
        drug2prot[r["drug_id"]].append(r["protein_id"])

    prot2gene = defaultdict(list)
    for r in pg:
        prot2gene[r["protein_id"]].append(r["gene_id"].upper())

    gene2hpo = defaultdict(list)
    for r in gh:
        gene2hpo[r["gene_id"].upper()].append(r["hpo_id"])

    actions = {}
    if ACTIONS.exists():
        for r in read_tsv(ACTIONS):
            actions[(r["drug_id"], r["protein_id"])] = r.get("action","")

    # HPO→PT + PT prior
    hpo2pt = {}
    for r in read_tsv(HPO2PT):
        hpo2pt[r["hpo_id"]] = r["meddra_pt"]
    pt_prior = { r["meddra_pt"]: float(r["prior"]) for r in read_csv(PTPRIOR) } if PTPRIOR.exists() else {}

    # CPIC drug→gene prior
    pgx = defaultdict(float)
    if PGX.exists():
        for r in read_csv(PGX):
            key = (r["drug_id"], r["gene_id"].upper())
            pgx[key] = max(pgx[key], float(r.get("cpic_weight", 0.0)))

    # score per drug per HPO
    for drug in sorted(set(r["drug_id"] for r in dt)):
        scores = defaultdict(float)

        # accumulate path weights
        for p in drug2prot[drug]:
            w_act = action_weight(actions.get((drug,p), ""))
            for g in prot2gene[p]:
                for h in gene2hpo[g]:
                    scores[h] += 1.0 * w_act

        # add CPIC contribution
        for (d,g), w in pgx.items():
            if d != drug: continue
            for h in gene2hpo.get(g,[]):
                scores[h] += BETA * w

        # add PT prior
        for h in list(scores.keys()):
            pt = hpo2pt.get(h)
            if pt:
                scores[h] += LAMBDA * pt_prior.get(pt, 0.0)

        # dump results
        out = OUTDIR/f"{drug}_ranked_hpo.csv"
        with out.open("w",encoding="utf-8",newline="") as fo:
            w = csv.writer(fo); w.writerow(["hpo_id","score"])
            for h, s in sorted(scores.items(), key=lambda kv: (-kv[1], kv[0])):
                w.writerow([h, f"{s:.6f}"])

        # small debug
        top3 = [h for h,_ in sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))[:3]]
        with (OUTDIR/f"debug_{drug}.txt").open("w",encoding="utf-8") as f:
            f.write("top3 HPO: " + ", ".join(top3) + "\n")

    print(f"Wrote ranked HPO lists to {OUTDIR}")

if __name__=="__main__":
    main()
