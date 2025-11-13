#!/usr/bin/env python3
# Build PT prior from SIDER meddra_all_se.tsv.gz (robust name/type detection; avoid CID/UMLS IDs)
import gzip, csv, math, re
from pathlib import Path

RAW = Path("data/raw/sider/meddra_all_se.tsv.gz")
OUT = Path("data/interim/mappings/pt_prior.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

if not RAW.exists():
    raise SystemExit(f"Missing {RAW}")

TYPES = {"PT","LLT","HLT","HLGT","SOC"}

def sniff_cols(path, max_scan=8000):
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        first = f.readline().rstrip("\n").split("\t")
        # Try header names first
        hdr = [c.strip().lower() for c in first]
        cand_name = {"name","concept_name","meddra_concept_name","side_effect_name"}
        cand_type = {"type","concept_type","meddra_concept_type","meddra_type"}
        name_idx = next((i for i,c in enumerate(hdr) if c in cand_name), None)
        type_idx = next((i for i,c in enumerate(hdr) if c in cand_type), None)
        if name_idx is not None and type_idx is not None:
            return name_idx, type_idx, True

        # Fallback: sniff by content
        rows = []
        for li, line in enumerate(f, start=1):
            r = line.rstrip("\n").split("\t")
            rows.append(r)
            if li >= max_scan: break
        if not rows:
            return None, None, False

        ncols = max(len(r) for r in rows)
        type_votes = [0]*ncols
        name_votes = [0]*ncols

        def is_type_token(x):
            return x.strip().upper() in TYPES

        def is_name_like(x):
            x = x.strip()
            if not re.search(r"[A-Za-z]", x): return False         # needs letters
            if re.match(r"^(CID|C)\d+$", x, re.I): return False     # avoid CID123.. / C123.. (UMLS)
            if len(x) < 3: return False
            return True

        for r in rows:
            for j in range(len(r)):
                v = r[j]
                if is_type_token(v): type_votes[j] += 1
                if is_name_like(v):
                    # weight: prefer multi-word / has lowercase letters
                    score = 1 + (1 if " " in v else 0) + (1 if re.search(r"[a-z]", v) else 0)
                    name_votes[j] += score

        type_idx = max(range(ncols), key=lambda j: type_votes[j])
        name_cands = [j for j in range(ncols) if j != type_idx]
        name_idx = max(name_cands, key=lambda j: name_votes[j]) if name_cands else None
        return name_idx, type_idx, False

name_idx, type_idx, used_header = sniff_cols(RAW)
if name_idx is None or type_idx is None:
    raise SystemExit("Could not detect name/type columns in SIDER.")

pt_count = {}
with gzip.open(RAW, "rt", encoding="utf-8", errors="ignore") as f:
    _ = f.readline()  # header or first row, already used for sniffing
    for line in f:
        r = line.rstrip("\n").split("\t")
        if type_idx >= len(r) or name_idx >= len(r):
            continue
        ctype = r[type_idx].strip().upper()
        name  = r[name_idx].strip()
        if ctype == "PT" and name:
            pt_count[name] = pt_count.get(name, 0) + 1

vals = [math.log1p(c) for c in pt_count.values()]
mx = max(vals) if vals else 1.0
rows = [{"meddra_pt": pt, "prior": round(math.log1p(cnt)/mx, 6)}
        for pt, cnt in pt_count.items()]
rows.sort(key=lambda x: -x["prior"])

with open(OUT, "w", encoding="utf-8", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=["meddra_pt","prior"])
    w.writeheader(); w.writerows(rows)

print(f"Wrote PT prior: {OUT} ({len(rows)} PTs) | used_header={used_header} name_idx={name_idx} type_idx={type_idx}")
