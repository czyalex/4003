from pyshacl import validate
from rdflib import Graph
data = Graph().parse("rdf/ae_kg.ttl", format="turtle")
sh   = Graph().parse("shacl/shapes.ttl", format="turtle")
conforms, results_graph, results_text = validate(data, shacl_graph=sh, inference="rdfs", abort_on_first=False)
print("Conforms:", conforms)
print(results_text)
