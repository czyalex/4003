#!/usr/bin/env python3
import csv, re
from pathlib import Path

ROOT = Path("."); M = ROOT/"data"/"interim"/"mappings"
MANUAL_DT = M/"drug_targets_manual.tsv"; MANUAL_PG = M/"protein_gene_manual.tsv"
DT = M/"drug_targets.csv"; PG = M/"protein_gene.csv"; GH = M/"gene_hpo.csv"

def read_tsv(p):
    # 兼容 UTF-8 with BOM 的 TSV
    with open(p, encoding="utf-8-sig") as fh:
        return [dict(zip(r.keys(), [v.strip() for v in r.values()]))
                for r in csv.DictReader(fh, delimiter="\t")]

def write_csv(p, rows, header):
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p,"w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f, fieldnames=header); w.writeheader(); w.writerows(rows)

def main():
    upat = re.compile(r"^[A-NR-Z][0-9]{5}$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$")

    # 1) drug_targets.csv
    dt_rows=[]
    for r in read_tsv(MANUAL_DT):
        drug=r.get("drug_id"); tgt=r.get("target_uniprot_or_name","")
        prot = tgt if upat.match(tgt) else tgt.upper()
        if drug and prot: dt_rows.append({"drug_id":drug,"protein_id":prot})
    write_csv(DT, dt_rows, ["drug_id","protein_id"])

    # 2) protein_gene.csv
    pg_rows=[]
    for r in read_tsv(MANUAL_PG):
        up=r.get("protein_uniprot"); gg=r.get("gene_symbol")
        if up and gg: pg_rows.append({"protein_id":up,"gene_id":gg})
    write_csv(PG, pg_rows, ["protein_id","gene_id"])

    # 3) gene_hpo.csv (初始化；后续用 HPO raw 扩充)
    gh_init=[{"gene_id":"HMGCR","hpo_id":"HP:0003198"},
             {"gene_id":"SCN1A","hpo_id":"HP:0012384"}]
    write_csv(GH, gh_init, ["gene_id","hpo_id"])

    print("Wrote:", DT, PG, GH)

if __name__=="__main__": main()
