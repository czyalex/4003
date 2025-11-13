#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fill missing UniProt -> gene_symbol mappings by querying UniProt REST API.
- Reads proteins from data/interim/mappings/drug_targets.csv
- Reads existing data/interim/mappings/protein_gene.csv
- Fills ONLY missing accessions
Outputs:
- Overwrites data/interim/mappings/protein_gene.csv with merged rows
- Writes a short report of unresolved accessions
"""
import csv, time, requests
from pathlib import Path

M = Path("data/interim/mappings")
DT = M/"drug_targets.csv"
PG = M/"protein_gene.csv"
REPORT = Path("reports")/("formal_"+__import__("datetime").date.today().isoformat())/"fill_missing_protein_gene.txt"

UNIPROT_ITEM = "https://rest.uniprot.org/uniprotkb/{}.json"
UNIPROT_TSV  = "https://rest.uniprot.org/uniprotkb/search?query=accession:{}&fields=gene_primary,gene_synonym&format=tsv"

def read_csv(p, encoding="utf-8-sig"):
    if not p.exists(): return []
    return [r for r in csv.DictReader(open(p, encoding=encoding))]

def write_csv(p, rows, header):
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8", newline="") as f:
        w=csv.DictWriter(f, fieldnames=header); w.writeheader(); w.writerows(rows)

def fetch_gene_from_item(acc):
    try:
        r = requests.get(UNIPROT_ITEM.format(acc), timeout=30)
        if r.status_code != 200: return ""
        j = r.json()
        genes = j.get("genes") or []
        if genes:
            g = genes[0]
            nm = (g.get("geneName") or {}).get("value") or ""
            if nm: return nm
            syns = g.get("synonyms") or []
            if syns:
                nm = (syns[0].get("value") or "")
                if nm: return nm
        return ""
    except Exception:
        return ""

def fetch_gene_from_tsv(acc):
    try:
        r = requests.get(UNIPROT_TSV.format(acc), timeout=30)
        if r.status_code != 200: return ""
        lines = r.text.strip().splitlines()
        if len(lines) >= 2:
            # Header: "Entry   Gene Names (primary)   Gene Names (synonym)"
            parts = lines[1].split("\t")
            if len(parts) >= 2 and parts[1].strip():
                return parts[1].strip()
            # fallback: synonym
            if len(parts) >= 3 and parts[2].strip():
                return parts[2].split(" ")[0].strip()
        return ""
    except Exception:
        return ""

def main():
    dt = read_csv(DT)
    pg = read_csv(PG)
    have = {(r["protein_id"] or "").strip(): (r["gene_id"] or "").strip() for r in pg if r.get("protein_id")}
    need = sorted({(r["protein_id"] or "").strip() for r in dt if r.get("protein_id")} - set(have.keys()))
    new_rows=[]
    unresolved=[]

    for i,acc in enumerate(need, 1):
        gene = fetch_gene_from_item(acc)
        if not gene:
            gene = fetch_gene_from_tsv(acc)
        if gene:
            new_rows.append({"protein_id": acc, "gene_id": gene})
        else:
            unresolved.append(acc)
        time.sleep(0.1)  # polite rate limit

    # merge & dedup
    merged = { (r["protein_id"], r["gene_id"]) for r in pg if r.get("protein_id") and r.get("gene_id") }
    for r in new_rows: merged.add((r["protein_id"], r["gene_id"]))
    out_rows = [{"protein_id":p, "gene_id":g} for (p,g) in sorted(merged)]

    write_csv(PG, out_rows, ["protein_id","gene_id"])

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write(f"Filled {len(new_rows)} accessions.\n")
        if unresolved:
            f.write("Unresolved:\n" + "\n".join(unresolved) + "\n")
    print(f"Updated {PG} ({len(out_rows)} rows). Unresolved: {len(unresolved)}")
    if unresolved:
        print("Unresolved list written to", REPORT)

if __name__ == "__main__":
    main()
