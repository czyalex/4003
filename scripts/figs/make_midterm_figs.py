# scripts/figs/make_midterm_figs.py
import os, re, csv, math, json, gzip
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(".")
FIG = ROOT/"figs"; FIG.mkdir(exist_ok=True, parents=True)

def latest_formal():
    dirs = sorted((ROOT/"reports").glob("formal_*"), key=lambda p: p.stat().st_mtime, reverse=True)
    return dirs[0] if dirs else None

def read_eval(formal_dir):
    rows = []
    with open(formal_dir/"eval_metrics.csv", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            rows.append({"drug": r["drug"], "P@10": float(r["P@10"]), "nDCG@10": float(r["nDCG@10"])})
    return rows

def plot_bar(rows, key, out):
    labs = [r["drug"] for r in rows]
    vals = [r[key] for r in rows]
    plt.figure()
    plt.bar(labs, vals)
    plt.ylabel(key)
    plt.title(key)
    plt.ylim(0, 1.0)
    for i,v in enumerate(vals):
        plt.text(i, v+0.02, f"{v:.3f}", ha="center")
    plt.savefig(out, bbox_inches="tight", dpi=160)
    plt.close()

def parse_hpo2pt_report():
    # 读取最新 hpo2pt_coverage.txt
    reps = sorted((ROOT/"reports").glob("formal_*/hpo2pt_coverage.txt"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not reps: return (0,0,0.0)
    t = reps[0].read_text(encoding="utf-8", errors="ignore")
    def get_num(pat):
        m = re.search(pat, t)
        return float(m.group(1)) if m else 0
    total = int(get_num(r"HPO nodes scanned:\s+(\d+)"))
    hits  = int(get_num(r"HPO->PT matches\s*:\s*(\d+)"))
    rate  = float(get_num(r"Match rate\s*:\s*([0-9.]+)"))
    return total, hits, rate

def plot_coverage(total, hits, rate, out):
    plt.figure()
    plt.bar(["Matched","Unmatched"], [hits, max(0,total-hits)])
    plt.title(f"HPO→PT coverage   matched={hits} / total={total}   rate={rate:.2%}")
    plt.ylabel("# HPO")
    plt.savefig(out, bbox_inches="tight", dpi=160)
    plt.close()

def targets_per_drug(out):
    dt = list(csv.DictReader(open(ROOT/"data/interim/mappings/drug_targets.csv", encoding="utf-8-sig")))
    from collections import Counter
    c = Counter(r["drug_id"] for r in dt)
    labs, vals = zip(*sorted(c.items()))
    plt.figure()
    plt.bar(labs, vals)
    plt.title("Targets per drug")
    plt.ylabel("# proteins")
    for i,v in enumerate(vals): plt.text(i, v+0.05, str(v), ha="center")
    plt.savefig(out, bbox_inches="tight", dpi=160)
    plt.close()

def main():
    formal = latest_formal()
    assert formal, "No reports/formal_* found."
    rows = read_eval(formal)
    plot_bar(rows, "P@10",   FIG/"p10_bar.png")
    plot_bar(rows, "nDCG@10",FIG/"ndcg_bar.png")
    total, hits, rate = parse_hpo2pt_report()
    plot_coverage(total, hits, rate, FIG/"hpo2pt_coverage.png")
    targets_per_drug(FIG/"targets_per_drug.png")
    print("Wrote figs to", FIG)

if __name__ == "__main__":
    main()
