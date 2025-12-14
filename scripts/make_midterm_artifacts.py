import os, csv
from pathlib import Path
from rdflib import Graph

report = Path("reports")/("midterm_" + __import__("datetime").date.today().isoformat())
report.mkdir(parents=True, exist_ok=True)

# 1) R1 路径导出为 CSV
g = Graph()
g.parse("rdf/schema.ttl", format="turtle")
g.parse("rdf/ae_kg.ttl",  format="turtle")
q = Path("queries/r1_paths.sparql").read_text(encoding="utf-8")
rows = [[str(x) for x in row] for row in g.query(q)]
with open(report/"r1_paths.csv","w",newline="",encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["drug","protein","gene","hpo"]); w.writerows(rows)

# 2) SHACL（尽力而为）
try:
    from pyshacl import validate
    sh = Graph().parse("shacl/shapes.ttl", format="turtle")
    conforms, _, results_text = validate(g, shacl_graph=sh, inference="rdfs", abort_on_first=False)
    (report/"shacl_report.txt").write_text(str(results_text), encoding="utf-8")
except Exception as e:
    (report/"shacl_report.txt").write_text("pyshacl not available or error: "+str(e), encoding="utf-8")

# 3) 三元组数 & 谓词使用频次
triples = list(g.triples((None,None,None)))
from collections import Counter
preds = Counter([str(p) for (_,p,_) in triples])
with open(report/"triple_counts.txt","w",encoding="utf-8") as f:
    f.write(f"Total triples: {len(triples)}\n")
with open(report/"predicate_usage.csv","w",newline="",encoding="utf-8") as f:
    w=csv.writer(f); w.writerow(["predicate","count"])
    for k,v in sorted(preds.items(), key=lambda x:-x[1]): w.writerow([k,v])

# 4) 评测（可选，存在就跑）
if Path("scripts/eval/run_eval.py").exists():
    os.system(f"python scripts/eval/run_eval.py > {report/'eval_stdout.txt'}")
print("Artifacts written to", report)
