#!/usr/bin/env python3
import gzip
from collections import Counter
from pathlib import Path

p = Path("data/raw/sider/meddra_all_se.tsv.gz")
assert p.exists(), "sider file not found: data/raw/sider/meddra_all_se.tsv.gz"

with gzip.open(p, "rt", encoding="utf-8", errors="ignore") as f:
    header = f.readline().rstrip("\n").split("\t")
    print("HEADER:", header)
    idx = {name.lower(): i for i,name in enumerate(header)}
    def pick(cands, fb=None):
        for c in cands:
            i = idx.get(c.lower())
            if i is not None: return i
        return fb
    i_type = pick(["concept_type","side_effect_type","meddra_type","type"], 4)
    i_name = pick(["concept_name","side_effect_name","meddra_name","name"], 3)
    i_drug = pick(["stitch_compound_id1","stitch_id1","drug_id","stitch_id"], 0)
    print("picked indices:", {"type":i_type, "name":i_name, "drug":i_drug})

    typ = Counter()
    examples = {}
    for n,line in enumerate(f, start=1):
        r = line.rstrip("\n").split("\t")
        t = (r[i_type] if i_type is not None and i_type < len(r) else "").strip()
        typ[t] += 1
        if t not in examples and i_name is not None and i_name < len(r):
            examples[t] = r[i_name]
        if n >= 200000: break

    print("concept_type counts (top 20):")
    for k,v in typ.most_common(20):
        print(f"  {repr(k)} : {v}")
    print("examples by type (sample):")
    for k in list(examples)[:10]:
        print(f"  {repr(k)} -> {examples[k]}")
