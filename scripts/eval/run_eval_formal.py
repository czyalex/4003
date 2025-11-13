#!/usr/bin/env python3
import csv, math, re
from pathlib import Path
from datetime import date

ROOT = Path(".")
GOLD = ROOT / "eval" / "gold"
MAP  = ROOT / "data" / "interim" / "mappings" / "hpo_meddra_map.tsv"

HP_PAT = re.compile(r"^HP:\d{7}$")

def read_txtlist(p):
    # 关键：utf-8-sig 读入 + 去掉显式 BOM
    lines = [line.strip() for line in open(p, encoding="utf-8-sig") if line.strip()]
    return [s.lstrip("\ufeff") for s in lines]

def p_at_k(ranked, gold, k=10):
    denom = max(1, min(k, len(ranked)))
    return sum(1 for x in ranked[:k] if x in gold) / float(denom)

def ndcg(ranked, gold, k=10):
    def _dcg(arr):
        s=0.0
        for i, x in enumerate(arr[:k], start=1):
            rel = 1.0 if x in gold else 0.0
            s += (2**rel - 1)/math.log2(i+1)
        return s
    ideal_hits = min(k, len(gold), max(1, len(ranked)))
    ideal = sum((2**1 - 1)/math.log2(i+1) for i in range(1, ideal_hits+1))
    return (_dcg(ranked)/ideal) if ideal>0 else 0.0

def load_map_tsv(p):
    m = {}
    if not p.exists(): return m
    with open(p, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            h = (r.get("hpo_id") or "").strip().upper()
            pt = (r.get("meddra_pt") or "").strip()
            if h and pt: m[h] = pt
    return m

def project_ranked_to_gold_space(ranked_hpo_ids, gold_items, hpo2pt, prefer_pt=True):
    gold_norm = [x.lstrip("\ufeff").strip() for x in gold_items]
    is_all_hpo = all(HP_PAT.match(x.upper()) for x in gold_norm)

    if prefer_pt and is_all_hpo:
        # 尝试把 gold 的 HPO→PT
        gold_pts = { (hpo2pt.get(x.upper()) or "").casefold() for x in gold_norm if hpo2pt.get(x.upper()) }
        ranked_pts = [ (hpo2pt.get(h.upper()) or "").casefold() for h in ranked_hpo_ids if hpo2pt.get(h.upper()) ]
        if gold_pts and ranked_pts:
            return ranked_pts, gold_pts, "PT"  # 覆盖到 PT 空间

    # 回退：HPO 空间
    if is_all_hpo:
        return [h.upper() for h in ranked_hpo_ids], {x.upper() for x in gold_norm}, "HPO"

    # gold 不是 HPO（认为是 PT 文本）
    ranked_pts = []
    for h in ranked_hpo_ids:
        pt = hpo2pt.get(h.upper())
        if pt: ranked_pts.append(pt.casefold())
    gold_pts = {x.casefold() for x in gold_norm}
    return ranked_pts, gold_pts, "PT"


def find_pred_dir():
    base = ROOT/"reports"
    cands = sorted(base.glob("formal_*"), key=lambda p: p.stat().st_mtime, reverse=True)
    for d in cands:
        if (d/"CHEMBL1064_ranked_hpo.csv").exists() or (d/"CHEMBL108_ranked_hpo.csv").exists():
            return d
    today = base/("formal_"+date.today().isoformat())
    today.mkdir(parents=True, exist_ok=True)
    return today

def main():
    out = find_pred_dir()
    out.mkdir(parents=True, exist_ok=True)
    drugs = [("CHEMBL1064","simvastatin"), ("CHEMBL108","carbamazepine")]

    hpo2pt = load_map_tsv(MAP)
    with open(out/"eval_debug.txt","w",encoding="utf-8") as d:
        d.write(f"Pred dir: {out}\n")
        d.write(f"Loaded HPO->PT mappings: {len(hpo2pt)}\n")
        for probe in ["HP:0003198","HP:0012384","HP:0000988","HP:0016266","HP:0005558"]:
            d.write(f"{probe} -> {hpo2pt.get(probe)}\n")

    with open(out/"eval_metrics.csv","w",encoding="utf-8",newline="") as fo:
        w = csv.writer(fo); w.writerow(["drug","k","P@10","nDCG@10","gold_space"])
        for chembl_id, disp in drugs:
            pred_path = out/f"{chembl_id}_ranked_hpo.csv"
            gold_path = GOLD/f"{disp}.txt"
            if not Path(pred_path).exists() or not gold_path.exists():
                continue
            ranked_hpo = [r["hpo_id"] for r in csv.DictReader(open(pred_path,encoding="utf-8-sig"))]
            gold_items = read_txtlist(gold_path)

            ranked_proj, gold_proj, space = project_ranked_to_gold_space(ranked_hpo, gold_items, hpo2pt)

            with open(out/f"debug_{disp}.txt","w",encoding="utf-8") as dd:
                dd.write("ranked (HPO) top3: "+", ".join(ranked_hpo[:3])+"\n")
                dd.write("ranked_proj (eval space) top3: "+", ".join(ranked_proj[:3])+"\n")
                dd.write("gold_proj sample: "+", ".join(list(gold_proj)[:5])+"\n")

            k=10
            w.writerow([disp, k, f"{p_at_k(ranked_proj,gold_proj,k):.3f}", f"{ndcg(ranked_proj,gold_proj,k):.3f}", space])
    print("Wrote eval_metrics.csv to", out)

if __name__=="__main__": main()
