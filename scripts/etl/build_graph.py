#!/usr/bin/env python3
"""Build a tiny RDF graph from minimal CSV mappings.

Inputs (place into data/interim/mappings/):
- drug_targets.csv: columns [drug_id, protein_id]
- protein_gene.csv:  columns [protein_id, gene_id]
- gene_hpo.csv:      columns [gene_id, hpo_id]

Output: rdf/ae_kg.ttl
"""
import pandas as pd
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, XSD

A = Namespace("http://example.org/ae-kg#")

def uri(kind, id_):
    return URIRef(f"http://example.org/{kind}/{id_}")

def main():
    g = Graph()
    g.bind("", A)
    # Load mappings
    dt = pd.read_csv("data/interim/mappings/drug_targets.csv")
    pg = pd.read_csv("data/interim/mappings/protein_gene.csv")
    gh = pd.read_csv("data/interim/mappings/gene_hpo.csv")
    # Types
    for _, r in dt.iterrows():
        g.add((uri("drug", r["drug_id"]),    RDF.type, A.Drug))
        g.add((uri("protein", r["protein_id"]), RDF.type, A.Protein))
        g.add((uri("drug", r["drug_id"]), A.actsOn, uri("protein", r["protein_id"])))
    for _, r in pg.iterrows():
        g.add((uri("protein", r["protein_id"]), RDF.type, A.Protein))
        g.add((uri("gene", r["gene_id"]), RDF.type, A.Gene))
        g.add((uri("protein", r["protein_id"]), A.encodedBy, uri("gene", r["gene_id"])))
    for _, r in gh.iterrows():
        g.add((uri("gene", r["gene_id"]), RDF.type, A.Gene))
        g.add((uri("phenotype", r["hpo_id"]), RDF.type, A.Phenotype))
        g.add((uri("gene", r["gene_id"]), A.hasPhenotype, uri("phenotype", r["hpo_id"])))
    g.serialize("rdf/ae_kg.ttl", format="turtle")
    print("Wrote rdf/ae_kg.ttl")
if __name__ == "__main__":
    main()
