#!/usr/bin/env python3
import gzip, re
from pathlib import Path

p = Path("data/raw/sider/meddra_all_se.tsv.gz")
assert p.exists(), f"not found: {p}"

uniq_pt=set(); uniq_drugs=set(); uniq_all=set()

with gzip.open(p, "rt", encoding="utf-8", errors="ignore") as f:
    first = f.readline().rstrip("\n").split("\t")
    # 判定是否“无表头”
    # 特征：某列就是 'PT' 或 'LLT'；第一列像 'CID\d+'
    headerless = any(x.strip().upper() in {"PT","LLT"} for x in first) or \
                 (len(first) >= 1 and re.match(r"^CID\d+$", first[0]) is not None)

    if headerless:
        # 你这份文件的结构：0=drug, 3=PT/LLT, 5=name（最后一列）
        def idx_type(tokens):
            for i,x in enumerate(tokens):
                if x.strip().upper() in {"PT","LLT"}: return i
            return 3
        i_drug = 0
        i_type = idx_type(first)
        i_name = len(first) - 1
        # 先处理首行
        cname = (first[i_name] if i_name < len(first) else "").strip()
        ctype = (first[i_type] if i_type < len(first) else "").strip().upper()
        drug  = (first[i_drug] if i_drug < len(first) else "").strip()
        if cname:
            uniq_all.add(cname)
            if ctype == "PT": uniq_pt.add(cname)
        if drug: uniq_drugs.add(drug)
        # 再处理剩余行
        for line in f:
            r = line.rstrip("\n").split("\t")
            cname = (r[i_name] if i_name < len(r) else "").strip()
            ctype = (r[i_type] if i_type < len(r) else "").strip().upper()
            drug  = (r[i_drug] if i_drug < len(r) else "").strip()
            if cname:
                uniq_all.add(cname)
                if ctype == "PT": uniq_pt.add(cname)
            if drug: uniq_drugs.add(drug)
    else:
        # 有表头（以防未来换文件），尽量兼容
        hdr = [h.lower() for h in first]
        idx = {name:i for i,name in enumerate(first)}
        def pick(cands, fb=None):
            for c in cands:
                j = idx.get(c)
                if j is not None: return j
            return fb
        i_type = pick(["concept_type","side_effect_type","meddra_type","type"], 4)
        i_name = pick(["concept_name","side_effect_name","meddra_name","name"], 3)
        i_drug = pick(["stitch_compound_id1","stitch_id1","drug_id","stitch_id"], 0)
        for line in f:
            r = line.rstrip("\n").split("\t")
            cname = (r[i_name] if i_name is not None and i_name < len(r) else "").strip()
            ctype = (r[i_type] if i_type is not None and i_type < len(r) else "").strip().upper()
            drug  = (r[i_drug] if i_drug is not None and i_drug < len(r) else "").strip()
            if cname:
                uniq_all.add(cname)
                if ctype == "PT": uniq_pt.add(cname)
            if drug: uniq_drugs.add(drug)

print("SIDER unique PT terms :", len(uniq_pt))
print("SIDER unique drugs    :", len(uniq_drugs))
print("SIDER ALL terms (any) :", len(uniq_all))
