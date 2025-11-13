#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fetch protein actions (mechanisms) from ChEMBL and resolve UniProt accessions.
Inputs:
  - drugs from CLI args or env CHEMBL_DRUGS="CHEMBL1064,CHEMBL108"
Outputs:
  - data/raw/chembl/chembl_mechanisms.csv                (raw enriched)
  - data/interim/mappings/drug_targets_enriched.csv      (drug_id, protein_id, action_type, is_moa, source)
  - data/interim/mappings/drug_targets.csv               (drug_id, protein_id)  # minimal for R1
"""
import os, sys, time, csv, requests
from pathlib import Path

ROOT = Path(".")
RAW  = ROOT/"data"/"raw"/"chembl"
OUTM = ROOT/"data"/"interim"/"mappings"
RAW.mkdir(parents=True, exist_ok=True)
OUTM.mkdir(parents=True, exist_ok=True)

MECH_API   = "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
TARGET_API = "https://www.ebi.ac.uk/chembl/api/data/target/{}.json"

def drugs_from_args_env():
    if len(sys.argv) > 1:
        return [x.strip() for x in sys.argv[1:] if x.strip()]
    env = os.environ.get("CHEMBL_DRUGS", "")
    if env:
        return [x.strip() for x in env.split(",") if x.strip()]
    # 默认用辛伐他汀+卡马西平
    return ["CHEMBL1064","CHEMBL108"]

def fetch_mechanisms(chembl_id, limit=200):
    out=[]
    offset=0
    while True:
        r = requests.get(MECH_API, params={"molecule_chembl_id": chembl_id, "limit": limit, "offset": offset}, timeout=60)
        r.raise_for_status()
        data = r.json()
        mechs = data.get("mechanisms") or data.get("mechanism") or []
        out.extend(mechs)
        if len(mechs) < limit: break
        offset += limit; time.sleep(0.1)
    return out

def resolve_accessions(target_chembl_id):
    """Fallback: query target endpoint to get target_components[].accession"""
    try:
        r = requests.get(TARGET_API.format(target_chembl_id), timeout=60)
        if r.status_code != 200:
            return []
        j = r.json()
        comps = j.get("target_components") or []
        accs=[]
        for c in comps:
            acc = (c.get("accession") or c.get("component_accession") or "").strip()
            if acc: accs.append(acc)
        return list(dict.fromkeys(accs))
    except Exception:
        return []

def flatten_rows(mechs):
    rows=[]
    for m in mechs:
        d = (m.get("molecule_chembl_id") or "").strip()
        tchembl = (m.get("target_chembl_id") or "").strip()
        action  = (m.get("action_type") or "").strip()
        moa_txt = (m.get("mechanism_of_action") or "").strip()
        is_moa  = 1 if moa_txt else 0

        # 先尝试机制里自带的 components
        comps = m.get("target_components") or []
        accs=[]
        for c in comps:
            acc = (c.get("accession") or "").strip()
            if acc: accs.append(acc)

        # 若没有 accession，回查 target 接口
        if not accs and tchembl:
            accs = resolve_accessions(tchembl)

        if accs:
            for acc in accs:
                rows.append({
                    "drug_id": d,
                    "protein_id": acc,
                    "target_chembl_id": tchembl,
                    "action_type": action,
                    "is_moa": is_moa,
                    "mechanism_of_action": moa_txt,
                    "source": "ChEMBL"
                })
        else:
            # 实在解析不到 accession，也保留一行（protein_id 留空，后续会被最小映射跳过）
            rows.append({
                "drug_id": d,
                "protein_id": "",
                "target_chembl_id": tchembl,
                "action_type": action,
                "is_moa": is_moa,
                "mechanism_of_action": moa_txt,
                "source": "ChEMBL"
            })
    return rows

def write_csv(path, rows, header):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader(); w.writerows(rows)

def main():
    drugs = drugs_from_args_env()
    all_rows=[]
    for d in drugs:
        mechs = fetch_mechanisms(d)
        rows  = flatten_rows(mechs)
        print(f"[{d}] mechanisms: {len(mechs)}  -> resolved rows: {sum(1 for r in rows if r['protein_id'])}")
        all_rows.extend(rows)
        time.sleep(0.2)

    raw_path = RAW/"chembl_mechanisms.csv"
    write_csv(raw_path, all_rows, ["drug_id","protein_id","target_chembl_id","action_type","is_moa","mechanism_of_action","source"])

    # 选“最好”的一条（优先 is_moa=1）
    best={}
    for r in all_rows:
        d=r["drug_id"]; p=r["protein_id"]
        if not p: continue
        key=(d,p)
        cur=best.get(key)
        if (cur is None) or (r["is_moa"]>cur["is_moa"]):
            best[key]=r

    enr_rows=[]
    for (d,p),r in best.items():
        enr_rows.append({"drug_id": d, "protein_id": p, "action_type": r["action_type"], "is_moa": r["is_moa"], "source": "ChEMBL"})

    write_csv(OUTM/"drug_targets_enriched.csv", enr_rows, ["drug_id","protein_id","action_type","is_moa","source"])
    write_csv(OUTM/"drug_targets.csv", [{"drug_id":r["drug_id"],"protein_id":r["protein_id"]} for r in enr_rows], ["drug_id","protein_id"])

    print("Wrote:", raw_path)
    print("Wrote:", OUTM/"drug_targets_enriched.csv")
    print("Wrote:", OUTM/"drug_targets.csv")

if __name__=="__main__":
    main()
