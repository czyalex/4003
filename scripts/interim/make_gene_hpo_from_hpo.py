#!/usr/bin/env python3
# 将 HPO 的 genes_to_phenotype.txt 转成 gene_hpo.csv
# 保留：结构路径基因（来自 protein_gene.csv）∪ CPIC 基因（来自 drug_gene_pgx.csv）

import csv
from pathlib import Path

ROOT = Path(".")
G2P = ROOT/"data/raw/hpo/genes_to_phenotype.txt"
PG  = ROOT/"data/interim/mappings/protein_gene.csv"
PGX = ROOT/"data/interim/mappings/drug_gene_pgx.csv"
OUT = ROOT/"data/interim/mappings/gene_hpo.csv"

def upper_set_from_csv_col(p: Path, col: str):
    s=set()
    if p.exists():
        with open(p, encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                v=(r.get(col) or "").strip().upper()
                if v: s.add(v)
    return s

def main():
    struct_genes = upper_set_from_csv_col(PG,  "gene_id")
    cpic_genes   = upper_set_from_csv_col(PGX, "gene_id")
    keep_genes   = struct_genes | cpic_genes

    pairs=set()
    with open(G2P, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            # entrez-id, gene-symbol, HPO-ID, HPO-name, ...
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3: continue
            gene = (parts[1] or "").strip().upper()
            hpo  = (parts[2] or "").strip().upper()
            if not hpo.startswith("HP:"): continue
            if keep_genes and gene not in keep_genes: continue
            pairs.add((gene, hpo))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["gene_id","hpo_id"])
        for g,h in sorted(pairs):
            w.writerow([g,h])

    print(f"Wrote {OUT} ({len(pairs)} rows). Included CPIC genes: {len(cpic_genes)}")

if __name__ == "__main__":
    main()
