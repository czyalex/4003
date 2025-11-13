#!/usr/bin/env python3
import csv, gzip, json, re
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(".")
OUTD = ROOT/"reports"/"inventory"
OUTD.mkdir(parents=True, exist_ok=True)

P = {
    "drug_targets": ROOT/"data/interim/mappings/drug_targets.csv",
    "protein_gene": ROOT/"data/interim/mappings/protein_gene.csv",
    "gene_hpo":    ROOT/"data/interim/mappings/gene_hpo.csv",
    "hpo_pt_map":  ROOT/"data/interim/mappings/hpo_meddra_map.tsv",
    "pgx":         ROOT/"data/interim/mappings/drug_gene_pgx.csv",
    "hp_json":     ROOT/"data/raw/hpo/hp.json",
    "hp_obo":      ROOT/"data/raw/hpo/hp.obo",
    "sider_se":    ROOT/"data/raw/sider/meddra_all_se.tsv.gz",
    "sider_drugs": ROOT/"data/raw/sider/drug_names.tsv.gz",
    "ttl":         ROOT/"rdf/ae_kg.ttl",
    "triple_counts": ROOT/"triple_counts.txt",
}

def read_csv_unique(path: Path, col: str, delimiter=","):
    if not path.exists(): return set()
    out=set()
    with open(path, encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f, delimiter=delimiter)
        for r in rdr:
            v=(r.get(col) or "").strip()
            if v: out.add(v)
    return out

def read_edges(path: Path, cols, delimiter=","):
    if not path.exists(): return 0
    cnt=0
    with open(path, encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f, delimiter=delimiter)
        for r in rdr:
            if all((r.get(c) or "").strip() for c in cols):
                cnt+=1
    return cnt

def count_hpo_total():
    """优先 hp.json；否则 hp.obo；若都失败，回退到 gene_hpo.csv 的去重"""
    # hp.json
    try:
        if P["hp_json"].exists():
            data = json.load(open(P["hp_json"], encoding="utf-8"))
            def walk(obj):
                if isinstance(obj, dict):
                    for k,v in obj.items():
                        if k=="id" and isinstance(v,str) and v.startswith("HP:"):
                            yield v
                        else:
                            yield from walk(v)
                elif isinstance(obj, list):
                    for x in obj: yield from walk(x)
            s=set(walk(data))
            if s: return len(s)
    except Exception:
        pass
    # hp.obo
    try:
        if P["hp_obo"].exists():
            n=0
            for line in open(P["hp_obo"], encoding="utf-8", errors="ignore"):
                if line.startswith("id: HP:"): n+=1
            if n>0: return n
    except Exception:
        pass
    # 兜底：gene_hpo.csv
    s = read_csv_unique(P["gene_hpo"], "hpo_id")
    return len(s)

def count_sider_pt():
    """统计 SIDER 的 PT 数量和涉及的 drug 数；自动适配 有表头/无表头 两种布局。"""
    import gzip
    se_path = P["sider_se"]
    if not se_path.exists():
        return 0, 0

    uniq_pt=set()
    uniq_drugs=set()

    with gzip.open(se_path, "rt", encoding="utf-8", errors="ignore") as f:
        first = f.readline().rstrip("\n").split("\t")

        # —— 识别是否为“无表头”版本：
        # 典型特征：len>=6 且 first[3] 恰好是 'PT' 或 'LLT'，first[0] 以 'CID' 开头
        headerless = (len(first) >= 6 and first[3] in ("PT", "LLT") and first[0].startswith("CID"))

        if headerless:
            # 无表头固定列位：0=drug, 3=type, 5=name
            i_drug, i_type, i_name = 0, 3, 5

            def handle_row(cols):
                if len(cols) <= max(i_drug, i_type, i_name): 
                    return
                drug = cols[i_drug].strip()
                typ  = cols[i_type].strip().upper()
                name = cols[i_name].strip()
                if drug: uniq_drugs.add(drug)
                if name and typ == "PT": uniq_pt.add(name)

            # 先处理“第一行数据”
            handle_row(first)
            # 再流式处理剩余各行
            for line in f:
                handle_row(line.rstrip("\n").split("\t"))

        else:
            # 有表头：按列名兼容多版本
            hdr = [h.strip() for h in first]
            idx = {name.lower(): i for i, name in enumerate(hdr)}
            def pick(cands, fb=None):
                for c in cands:
                    i = idx.get(c.lower())
                    if i is not None: return i
                return fb
            i_type = pick(["concept_type","side_effect_type","meddra_type","type"], 4)
            i_name = pick(["concept_name","side_effect_name","meddra_name","name"], 3)
            i_drug = pick(["stitch_compound_id1","stitch_id1","drug_id","stitch_id"], 0)

            def is_pt(s: str) -> bool:
                if not s: return False
                u = s.strip().upper()
                return (u == "PT") or (u.startswith("PT")) or ("PREFERRED" in u and "TERM" in u)

            for line in f:
                r = line.rstrip("\n").split("\t")
                ctype = (r[i_type] if i_type is not None and i_type < len(r) else "")
                cname = (r[i_name] if i_name is not None and i_name < len(r) else "")
                drug  = (r[i_drug] if i_drug is not None and i_drug < len(r) else "")
                if drug: uniq_drugs.add(drug.strip())
                if cname and is_pt(ctype): uniq_pt.add(cname.strip())

    return len(uniq_pt), len(uniq_drugs)


# — 实体规模 —
drug_ids     = read_csv_unique(P["drug_targets"], "drug_id")
proteins     = read_csv_unique(P["drug_targets"], "protein_id")
genes_pg     = {g.upper() for g in read_csv_unique(P["protein_gene"], "gene_id")}
genes_pgx    = {g.upper() for g in read_csv_unique(P["pgx"], "gene_id")}
genes_total  = genes_pg | genes_pgx
hpo_total    = count_hpo_total()
pt_total, sider_drug_total = count_sider_pt()
triples      = count_triples()

# — 边 —
edges_dt  = read_edges(P["drug_targets"], ["drug_id","protein_id"])
edges_pg  = read_edges(P["protein_gene"], ["protein_id","gene_id"])
edges_gh  = read_edges(P["gene_hpo"],    ["gene_id","hpo_id"])
edges_hppt= read_edges(P["hpo_pt_map"],  ["hpo_id","meddra_pt"], delimiter="\t")

pt_total, sider_drug_total, se_all_total, pt_lower = count_sider_pt()

inv_rows = [
    {"entity":"Drugs (MVP)",      "count":len(drug_ids),    "source":"You selected (CHEMBL IDs)",      "files":str(P["drug_targets"].relative_to(ROOT))},
    {"entity":"Proteins",         "count":len(proteins),    "source":"UniProt (via ChEMBL/DrugCentral)","files":str(P["drug_targets"].relative_to(ROOT))},
    {"entity":"Genes (struct)",   "count":len(genes_pg),    "source":"HGNC (via UniProt mapping)",      "files":str(P["protein_gene"].relative_to(ROOT))},
    {"entity":"Genes (CPIC)",     "count":len(genes_pgx),   "source":"CPIC prior",                      "files":str(P["pgx"].relative_to(ROOT))},
    {"entity":"Genes (union)",    "count":len(genes_total), "source":"struct + CPIC",                   "files":"protein_gene.csv + drug_gene_pgx.csv"},
    {"entity":"HPO terms (total)","count":hpo_total,        "source":"HPO hp.json / hp.obo / fallback", "files":"data/raw/hpo/hp.json | hp.obo | gene_hpo.csv"},
    {"entity":"MedDRA PT (SIDER)","count":pt_total,         "source":"SIDER (meddra_all_se.tsv.gz)",    "files":str(P["sider_se"].relative_to(ROOT)) if P["sider_se"].exists() else "N/A"},
]
# 若严格 PT 仍为 0，则附加两行“回退统计”避免 PPT 空缺
if pt_total == 0:
    inv_rows += [
        {"entity":"SIDER terms (all types)","count":se_all_total, "source":"SIDER (fallback all names)","files":str(P["sider_se"].relative_to(ROOT))},
        {"entity":"PT (≥ HPO↔SIDER lower bound)","count":pt_lower,"source":"Intersection with HPO→PT","files":"hpo_meddra_map.tsv ∩ SIDER"},
    ]

with open(OUTD/"data_inventory.csv","w",encoding="utf-8",newline="") as fo:
    w=csv.DictWriter(fo, fieldnames=["entity","count","source","files"])
    w.writeheader(); w.writerows(inv_rows)

with open(OUTD/"edge_counts.csv","w",encoding="utf-8",newline="") as fo:
    w=csv.DictWriter(fo, fieldnames=["mapping","edges","files"])
    w.writeheader()
    w.writerow({"mapping":"Drug→Protein","edges":edges_dt,"files":str(P["drug_targets"].relative_to(ROOT))})
    w.writerow({"mapping":"Protein→Gene","edges":edges_pg,"files":str(P["protein_gene"].relative_to(ROOT))})
    w.writerow({"mapping":"Gene→HPO","edges":edges_gh,"files":str(P["gene_hpo"].relative_to(ROOT))})
    w.writerow({"mapping":"HPO→PT","edges":edges_hppt,"files":str(P["hpo_pt_map"].relative_to(ROOT))})

with open(OUTD/"SUMMARY.txt","w",encoding="utf-8") as f:
    f.write("Data inventory summary\n")
    f.write(f"- Drugs (MVP): {len(drug_ids)}\n")
    f.write(f"- Proteins: {len(proteins)}\n")
    f.write(f"- Genes (struct ∪ CPIC): {len(genes_total)}\n")
    f.write(f"- HPO terms (total): {hpo_total}\n")
    f.write(f"- MedDRA PT (SIDER): {pt_total}\n")
    f.write(f"- RDF triples: {triples}\n")

# 简单柱状图
labels = [r["entity"] for r in inv_rows]
values = [r["count"]  for r in inv_rows]
plt.figure()
plt.bar(labels, values)
plt.xticks(rotation=45, ha="right")
plt.ylabel("Count")
plt.title("Project data inventory")
plt.tight_layout()
plt.savefig(OUTD/"data_inventory.png", dpi=160)
plt.close()

print("Wrote:", OUTD/"data_inventory.csv", OUTD/"edge_counts.csv", OUTD/"data_inventory.png", OUTD/"SUMMARY.txt")
