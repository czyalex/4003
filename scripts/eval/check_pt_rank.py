#!/usr/bin/env python3
import csv
from pathlib import Path

ROOT = Path(".")

def load_hpo2pt(path):
    """Load HPO -> MedDRA PT mapping from TSV."""
    mapping = {}
    with open(path, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            h = (r.get("hpo_id") or "").strip().upper()
            # be robust to column naming: pt_name or meddra_pt
            pt = (r.get("pt_name") or r.get("meddra_pt") or "").strip()
            if h and pt:
                mapping[h] = pt.casefold()
    return mapping

def pt_rank(pred_path, target_pts, hpo2pt):
    """Print the rank of given PT names in a ranked HPO prediction file."""
    ranked_hpo = []
    with open(pred_path, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for r in reader:
            ranked_hpo.append((r.get("hpo_id") or "").strip().upper())

    # project HPO sequence to PT sequence using the mapping
    ranked_pts = [hpo2pt.get(h) for h in ranked_hpo if hpo2pt.get(h)]

    # first occurrence rank for each PT
    seen = {}
    for i, pt in enumerate(ranked_pts, start=1):
        if pt not in seen:
            seen[pt] = i

    print(f"\n=== {Path(pred_path)} ===")
    for t in target_pts:
        t_norm = t.casefold()
        print(f"{t}: rank {seen.get(t_norm)}")

def main():
    hpo2pt = load_hpo2pt(ROOT / "data" / "interim" / "mappings" / "hpo_meddra_map.tsv")

    # ----- No CPIC: use the older run (formal_2025-11-12) -----
    pt_rank(
        ROOT / "reports" / "formal_2025-11-12" / "CHEMBL108_ranked_hpo.csv",
        ["Stevens-Johnson syndrome", "Toxic epidermal necrolysis"],
        hpo2pt,
    )
    pt_rank(
        ROOT / "reports" / "formal_2025-11-12" / "CHEMBL1064_ranked_hpo.csv",
        ["Myopathy"],
        hpo2pt,
    )

    # ----- With CPIC: use the new run (formal_2025-12-01) -----
    pt_rank(
        ROOT / "reports" / "formal_2025-12-01" / "CHEMBL108_ranked_hpo.csv",
        ["Stevens-Johnson syndrome", "Toxic epidermal necrolysis"],
        hpo2pt,
    )
    pt_rank(
        ROOT / "reports" / "formal_2025-12-01" / "CHEMBL1064_ranked_hpo.csv",
        ["Myopathy"],
        hpo2pt,
    )

if __name__ == "__main__":
    main()
