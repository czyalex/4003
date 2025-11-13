# scripts/interim/apply_gene_hpo_overrides.py
# Purpose: apply manual Gene->HPO overrides from a TSV with columns:
#   gene_symbol<TAB>hpo
# where "hpo" can be an HPO ID (e.g., HP:0003326) OR a human-readable
# HPO/MedDRA PT name (e.g., "Myalgia", "Stevens-Johnson syndrome").
# The script appends deduped overrides into data/interim/mappings/gene_hpo.csv.

import csv, json, os, sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
MAPP = BASE / "data" / "interim" / "mappings"
RAWHPO = BASE / "data" / "raw" / "hpo"
OVR = MAPP / "gene_hpo_overrides.tsv"
OUT = MAPP / "gene_hpo.csv"
HPO2PT = MAPP / "hpo_meddra_map.tsv"  # for resolving PT names -> HPO IDs
HP_JSON = RAWHPO / "hp.json"          # optional for resolving HPO labels -> IDs

def read_tsv(path):
    # Use utf-8-sig to swallow BOM if present
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

def load_gene_hpo():
    if not OUT.exists():
        return []
    with open(OUT, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))

def load_hpo_label_index():
    """Return dict: lowercase HPO label -> HPO ID (from hp.json if available)."""
    idx = {}
    if HP_JSON.exists():
        data = json.load(open(HP_JSON, "r", encoding="utf-8"))
        for node in data.get("graphs",[{}])[0].get("nodes",[]):
            if node.get("id","").startswith("HP:"):
                lbl = node.get("lbl")
                if lbl:
                    idx[lbl.lower()] = node["id"]
                # include exact synonyms if present
                for syn in node.get("meta",{}).get("synonyms",[]):
                    s = syn.get("val")
                    if s:
                        idx[s.lower()] = node["id"]
    return idx

def load_pt_to_hpo():
    """Return dict: lowercase MedDRA PT -> HPO ID (from hpo_meddra_map.tsv)."""
    d = {}
    if HPO2PT.exists():
        with open(HPO2PT, "r", encoding="utf-8", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                hpo = r.get("hpo_id","").strip()
                pt  = r.get("meddra_pt","").strip()
                if hpo.startswith("HP:") and pt:
                    d.setdefault(pt.lower(), set()).add(hpo)
    return d

def normalize_hpo(term, hpo_label_idx, pt2hpo):
    """Return an HPO ID for a given term."""
    t = term.strip()
    if not t:
        return None
    # If already an HPO ID
    if t.upper().startswith("HP:"):
        return t.upper()
    # Try HPO label
    hid = hpo_label_idx.get(t.lower())
    if hid:
        return hid
    # Try MedDRA PT -> HPO (from our derived map)
    cands = pt2hpo.get(t.lower())
    if cands:
        # pick the first deterministic ID
        return sorted(cands)[0]
    return None

def main():
    if not OVR.exists():
        print(f"Overrides file not found: {OVR}")
        sys.exit(1)

    overrides = read_tsv(OVR)
    # Defensive: cope with BOMed headers by probing keys
    def vget(row, *keys):
        for k in keys:
            if k in row: return row[k]
        # try stripping BOM from keys
        for k in row:
            if k.lstrip("\ufeff") in keys:
                return row[k]
        return None

    base = load_gene_hpo()
    base_pairs = {(r["gene_id"].upper().strip(), r["hpo_id"].upper().strip())
                  for r in base if r.get("gene_id") and r.get("hpo_id")}

    hpo_label_idx = load_hpo_label_index()
    pt2hpo = load_pt_to_hpo()

    added = 0
    for r in overrides:
        g = (vget(r, "gene_symbol") or "").upper().strip()
        term = (vget(r, "hpo") or "").strip()
        if not g or not term:
            continue
        hpo_id = normalize_hpo(term, hpo_label_idx, pt2hpo)
        if not hpo_id:
            print(f"[WARN] Cannot resolve '{term}' to an HPO ID; skipping")
            continue
        pair = (g, hpo_id)
        if pair not in base_pairs:
            base.append({"gene_id": g, "hpo_id": hpo_id})
            base_pairs.add(pair)
            added += 1

    write_csv(OUT, base, fieldnames=["gene_id","hpo_id"])
    print(f"Applied overrides: +{added} pairs. Total rows now: {len(base)}")
    print(f"Wrote {OUT}")

if __name__ == "__main__":
    main()
