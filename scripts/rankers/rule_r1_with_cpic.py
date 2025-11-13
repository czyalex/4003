#!/usr/bin/env python3
# R1 with CPIC + multi-path support + lambda * PT prior
import csv
from pathlib import Path
from collections import defaultdict

BETA   = 0.6   # CPIC 鏉冮噸锛堝悗缁彲璋冨弬锛?
LAMBDA = 0.3   # PT 鍏堥獙鏉冮噸锛堝悗缁彲璋冨弬锛?
K      = 10

ROOT = Path(".")
M    = ROOT/"data"/"interim"/"mappings"
RPT  = ROOT/"reports"/f"formal_{__import__('datetime').datetime.now():%Y-%m-%d}"
RPT.mkdir(parents=True, exist_ok=True)

def read_csv(path, delim=","):
    with open(path, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter=delim))

def write_rank(drug_id, hpo2score):
    rows = sorted(hpo2score.items(), key=lambda x:-x[1])
    with open(RPT/f"{drug_id}_ranked_hpo.csv","w",encoding="utf-8",newline="") as fo:
        w=csv.writer(fo); w.writerow(["hpo_id","score"]); w.writerows(rows)

def main():
    dt  = read_csv(M/"drug_targets.csv")            # drug_id, protein_id
    pg  = read_csv(M/"protein_gene.csv")            # protein_id, gene_id
    gh  = read_csv(M/"gene_hpo.csv")                # gene_id, hpo_id
    pgx = read_csv(M/"drug_gene_pgx.csv")           # drug_id, gene_id, cpic_weight
    hm  = read_csv(M/"hpo_meddra_map.tsv", "\t")    # hpo_id, meddra_pt

    # PT prior
    pt_prior = {r["meddra_pt"].lower(): float(r["prior"]) for r in read_csv(M/"pt_prior.csv")} if (M/"pt_prior.csv").exists() else {}

    # indices
    drug2prot = defaultdict(set)
    for r in dt: drug2prot[r["drug_id"]].add(r["protein_id"])

    prot2gene = defaultdict(set)
    for r in pg: prot2gene[r["protein_id"]].add(r["gene_id"].upper())

    gene2hpo  = defaultdict(set)
    for r in gh: gene2hpo[r["gene_id"].upper()].add(r["hpo_id"])

    hpo2pt    = {}
    for r in hm:
        if r.get("hpo_id") and r.get("meddra_pt"):
            hpo2pt[r["hpo_id"]] = r["meddra_pt"]

    drug2pgxgenes = defaultdict(list)
    for r in pgx:
        try:
            drug2pgxgenes[r["drug_id"]].append((r["gene_id"].upper(), float(r.get("cpic_weight",1.0))))
        except: pass

    # rank per drug
    for drug in sorted({r["drug_id"] for r in dt}):
        # collect reachable genes via structure
        genes_struct = set()
        for p in drug2prot[drug]:
            genes_struct |= prot2gene[p]

        # accumulate supports per HPO锛氱敱澶氬皯 unique genes 鏀寔
        hpo_supports = defaultdict(int)
        for g in genes_struct:
            for h in gene2hpo.get(g, []):
                hpo_supports[h] += 1

        # 褰掍竴鍖栵紙閬垮厤鍏?1.0锛?
        max_sup = max(hpo_supports.values()) if hpo_supports else 1
        hpo_score = {h: s/max_sup for h,s in hpo_supports.items()}

        # CPIC 鍙犲姞锛堜笉鏄?max锛?
        for g, w in drug2pgxgenes.get(drug, []):
            for h in gene2hpo.get(g, []):
                hpo_score[h] = hpo_score.get(h, 0.0) + BETA * w

        # PT prior 鍙犲姞
        for h in list(hpo_score.keys()):
            pt = hpo2pt.get(h, "")
            if pt:
                hpo_score[h] += LAMBDA * pt_prior.get(pt.lower(), 0.0)

        write_rank(drug, hpo_score)

    print(f"Wrote ranked HPO lists to {RPT}")

if __name__ == "__main__":
    main()

