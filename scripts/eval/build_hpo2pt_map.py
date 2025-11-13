#!/usr/bin/env python3
import json, csv, gzip, re, sys
from pathlib import Path
from datetime import date

ROOT = Path(".")
HPO_DIR   = ROOT/"data"/"raw"/"hpo"
SIDER_DIR = ROOT/"data"/"raw"/"sider"
OUT_MAP = ROOT/"data"/"interim"/"mappings"/"hpo_meddra_map.tsv"
OUT_RPT = ROOT/"reports"/("formal_"+date.today().isoformat())/"hpo2pt_coverage.txt"

_norm_rx = re.compile(r"[^a-z0-9 ]+")
def norm(s: str) -> str:
    s = s.casefold()
    s = _norm_rx.sub(" ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s

def open_text(p: Path):
    if str(p).endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8")
    return open(p, "r", encoding="utf-8")

def load_pt_vocab():
    terms = set()
    for fn in ("meddra_all_se.tsv", "meddra_all_se.tsv.gz"):
        p = SIDER_DIR / fn
        if not p.exists(): continue
        with open_text(p) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols:
                    pt = cols[-1].strip()
                    if pt: terms.add(pt)
    return terms

def parse_hp_json(p: Path):
    """return dict[HPO_ID] -> set(normalized strings) from hp.json (OBO Graph JSON)"""
    try:
        obj = json.load(open(p, "r", encoding="utf-8"))
    except Exception:
        return {}
    graphs = obj.get("graphs") or []
    if not graphs: 
        return {}
    nodes = graphs[0].get("nodes") or []
    hpo2strings = {}
    for n in nodes:
        nid = n.get("id","")
        if not nid.startswith("HP:"): 
            continue
        # skip property nodes
        if n.get("type") == "PROPERTY":
            continue
        strings = set()
        if n.get("lbl"): strings.add(n["lbl"])
        syns = (((n.get("meta") or {}).get("synonyms")) or [])
        for s in syns:
            v = (s.get("val") or "").strip()
            if v: strings.add(v)
        if strings:
            hpo2strings[nid] = {norm(x) for x in strings}
    return hpo2strings

def parse_hp_obo(p: Path):
    """fallback: parse hp.obo to labels/synonyms"""
    hpo2strings = {}
    if not p.exists(): return hpo2strings
    with open_text(p) as f:
        cur = None
        for line in f:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = {"id":None, "name":None, "syn":set()}
            elif not line and cur and cur["id"]:
                strings = set()
                if cur["name"]: strings.add(cur["name"])
                strings |= cur["syn"]
                if strings:
                    hpo2strings[cur["id"]] = {norm(x) for x in strings}
                cur = None
            elif cur is not None:
                if line.startswith("id: HP:"):
                    cur["id"] = line.split("id: ")[1].strip()
                elif line.startswith("name: "):
                    cur["name"] = line.split("name: ",1)[1].strip()
                elif line.startswith("synonym: "):
                    # synonym: "Some text" EXACT [...]
                    m = re.search(r'^synonym:\s+"(.+?)"\s', line)
                    if m: cur["syn"].add(m.group(1))
        # flush last term if file doesn't end with blank line
        if cur and cur["id"]:
            strings = set()
            if cur["name"]: strings.add(cur["name"])
            strings |= cur["syn"]
            if strings:
                hpo2strings[cur["id"]] = {norm(x) for x in strings}
    return hpo2strings

def main():
    # 1) load SIDER PT vocab
    pt_vocab = load_pt_vocab()
    if not pt_vocab:
        raise FileNotFoundError("No SIDER PT terms found under data/raw/sider/")
    pt_norm = {norm(x): x for x in pt_vocab}

    # 2) load HPO labels/synonyms (json first, then obo fallback)
    hpo2strings = {}
    hp_json = HPO_DIR/"hp.json"
    if hp_json.exists():
        hpo2strings = parse_hp_json(hp_json)
    if not hpo2strings:
        hp_obo = HPO_DIR/"hp.obo"
        if hp_obo.exists():
            hpo2strings = parse_hp_obo(hp_obo)

    total_hpo = len(hpo2strings)
    hits = []
    if total_hpo:
        for hpo_id, strings in hpo2strings.items():
            for s in strings:
                if s in pt_norm:
                    hits.append((hpo_id, pt_norm[s]))
                    break  # pick the first match

    OUT_MAP.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_MAP, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["hpo_id","meddra_pt"])
        for h, pt in hits:
            w.writerow([h, pt])

    OUT_RPT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_RPT, "w", encoding="utf-8") as r:
        r.write(f"HPO nodes scanned: {total_hpo}\n")
        r.write(f"HPO->PT matches : {len(hits)}\n")
        rate = (len(hits)/total_hpo) if total_hpo else 0.0
        r.write(f"Match rate      : {rate:.4f}\n")
        r.write(f"Output map      : {OUT_MAP}\n")

    print("Wrote:", OUT_MAP, "and", OUT_RPT)

if __name__ == "__main__":
    main()
