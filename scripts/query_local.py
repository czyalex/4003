from rdflib import Graph
qfile = "queries/r1_paths.sparql"
g = Graph()
g.parse("rdf/schema.ttl", format="turtle")
g.parse("rdf/ae_kg.ttl",  format="turtle")
with open(qfile, "r", encoding="utf-8") as f:
    q = f.read()
print("Running:", qfile)
for row in g.query(q):
    print([str(x) for x in row])
