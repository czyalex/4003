#!/usr/bin/env python3
"""Download small raw snapshots for MVP.

This is a stub with example endpoints to fill.
Usage:
  python scripts/download_raw.py --drug simvastatin --drug carbamazepine
"""
import argparse, pathlib, time, json, sys
from pathlib import Path

def main(drugs):
    root = Path("data/raw")
    for src in ["hpo","hgnc","uniprot","chembl","drugcentral","sider","cpic"]:
        for d in drugs:
            p = root / src / d.lower()
            p.mkdir(parents=True, exist_ok=True)
    # TODO: add real downloads (HPO, HGNC, etc.). For now we just drop README placeholders.
    for d in drugs:
        for src in ["hpo","hgnc","uniprot","chembl","drugcentral","sider","cpic"]:
            (root/src/d.lower()/"README.txt").write_text(f"Placeholder for {src}/{d}\n", encoding="utf-8")
    print("Created placeholder folders under data/raw/. Fill actual downloads next.")
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", action="append", required=True, help="Drug name")
    args = ap.parse_args()
    main(args.drug)
