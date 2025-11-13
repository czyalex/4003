#!/usr/bin/env python3
# Keep only HPO terms under "Phenotypic abnormality" (HP:0000118).
# Optionally keep leaf-only terms to avoid very broad parents.

import json, csv
from pathlib import Path
from collections import defaultdict, deque

BASE = Path("data/interim/mappings")
INP  = BASE/"gene_hpo.csv"
OUT  = BASE/"gene_hpo.filtered.csv"
HPJSON = Path("data/raw/hpo/hp.json")

KEEP_ROOT = "HP:0000118"  # Phenotypic abnormality
LEAF_ONLY = True          # set False if you want parents too (no discount in this script)

def build_graph():
    obj = json.loads(HPJSON.read_text(encoding="utf-8"))
    nodes = obj["graphs"][0]["nodes"]; edges = obj["graphs"][0]["edges"]
    parents = defaultdict(list); children = defaultdict(list)
    for e in edges:
        s = e.get("sub","").replace("_",":"); o = e.get("obj","").replace("_",":")
        if not s.startswith("HP:") or not o.startswith("HP:"): continue
        # is_a only
        if e.get("pred","").endswith("is_a"):
            parents[s].append(o); children[o].append(s)
    return parents, children

def is_descendant_of(parents, node, root):
    # BFS to root
    q = deque([node]); seen=set()
    while q:
        x = q.popleft()
        if x in seen: continue
        seen.add(x)
        if x == root: return True
        for p in parents.get(x,[]):
            q.append(p)
    return False

def main():
    parents, children = build_graph()
    rows = list(csv.DictReader(INP.open("r",encoding="utf-8")))
    keep = []
    # precompute leafs
    non_leaf = set()
    for _, kids in children.items():
        for k in kids: non_leaf.add(k)

    for r in rows:
        h = r["hpo_id"]
        if not is_descendant_of(parents, h, KEEP_ROOT):
            continue
        if LEAF_ONLY and h in non_leaf:
            continue
        keep.append(r)

    with OUT.open("w",encoding="utf-8",newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=["gene_id","hpo_id"])
        w.writeheader(); w.writerows(keep)

    print(f"Filtered {len(rows)} -> {len(keep)} rows. Output: {OUT}")

if __name__=="__main__":
    main()
