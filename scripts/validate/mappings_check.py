#!/usr/bin/env python3
import re, csv
from pathlib import Path

M = Path("data/interim/mappings")
problems=[]

def check_csv(p, required_cols):
    if not p.exists():
        problems.append(("MISSING_FILE", str(p)))
        return []
    # 兼容 UTF-8 with BOM，并在空文件时给出提示而不崩
    rows = list(csv.DictReader(open(p, encoding="utf-8-sig")))
    if not rows:
        problems.append(("EMPTY_FILE", str(p)))
        return []
    for col in required_cols:
        if col not in rows[0]:
            problems.append(("MISSING_COL", f"{p}::{col}"))
    return rows

def main():
    upat = re.compile(r"^[A-NR-Z][0-9]{5}$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$")

    dt = M/"drug_targets.csv"
    rows = check_csv(dt, ["drug_id","protein_id"])
    for i,r in enumerate(rows):
        if not upat.match(r["protein_id"]):
            problems.append(("BAD_UNIPROT", f"{dt} line {i+2} :: {r['protein_id']}"))

    pg = M/"protein_gene.csv"
    rows = check_csv(pg, ["protein_id","gene_id"])
    for i,r in enumerate(rows):
        if not upat.match(r["protein_id"]):
            problems.append(("BAD_UNIPROT", f"{pg} line {i+2} :: {r['protein_id']}"))

    gh = M/"gene_hpo.csv"
    rows = check_csv(gh, ["gene_id","hpo_id"])
    for i,r in enumerate(rows):
        # 修正后的 HPO 正则
        if not re.match(r"^HP:\d{7}$", r["hpo_id"]):
            problems.append(("BAD_HPO", f"{gh} line {i+2} :: {r['hpo_id']}"))

    pgx = M/"drug_gene_pgx.csv"
    if pgx.exists():
        with open(pgx, encoding="utf-8-sig") as f:
            rows=list(csv.DictReader(f))
        for i,r in enumerate(rows):
            try:
                float(r.get("cpic_weight","0") or 0.0)
            except:
                problems.append(("BAD_WEIGHT", f"{pgx} line {i+2} :: {r.get('cpic_weight')}"))

    out = Path("reports")/("formal_"+__import__("datetime").date.today().isoformat())
    out.mkdir(parents=True, exist_ok=True)
    op = out/"mapping_validation.txt"
    if problems:
        op.write_text("\n".join([f"{k}\t{v}" for k,v in problems]), encoding="utf-8")
        print("Validation: FAIL — see", op)
    else:
        op.write_text("OK", encoding="utf-8")
        print("Validation: OK")

if __name__=="__main__": main()
