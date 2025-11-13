#!/usr/bin/env python3
import csv, gzip
from pathlib import Path
from collections import defaultdict

RAW=Path("data/raw/sider"); OUT=Path("eval/gold")

def open_any(p): 
    return gzip.open(p, "rt", encoding="utf-8") if str(p).endswith(".gz") else open(p,"r",encoding="utf-8")

def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # 名称→STITCH ID
    name_to_ids=defaultdict(set)
    for fn in ["drug_names.tsv","drug_names.tsv.gz"]:
        p=RAW/fn
        if not p.exists(): continue
        for row in csv.reader(open_any(p), delimiter="\t"):
            if len(row)>=2: name_to_ids[row[1].lower()].add(row[0])

    targets=["simvastatin","carbamazepine"]
    target_ids=set().union(*(name_to_ids.get(t,set()) for t in targets))

    # 收集 AE 词
    drug_to_terms={t:set() for t in targets}
    for fn in ["meddra_all_se.tsv","meddra_all_se.tsv.gz"]:
        p=RAW/fn
        if not p.exists(): continue
        for row in csv.reader(open_any(p), delimiter="\t"):
            if len(row)<2: continue
            sid=row[0]; term=row[-1].strip()
            if sid not in target_ids: continue
            for d in targets:
                if sid in name_to_ids.get(d,set()): drug_to_terms[d].add(term)

    # 写 gold
    for d in targets:
        with open(OUT/f"{d}.txt","w",encoding="utf-8") as f:
            for t in sorted(drug_to_terms[d]): f.write(t+"\n")
    print("Gold written to", OUT)

if __name__=="__main__": main()
