# Midterm Checkpoint — AE-KG MVP (Toy Mapping)

## 1. Problem & Goal
- Pain: AE retrieval is stringy & opaque
- Goal: Auditable path rules + CPIC priors improve top-k AE retrieval

## 2. Data & Ontology
- Sources: HPO, HGNC, UniProt, ChEMBL/DrugCentral, SIDER, CPIC
- Minimal schema: Drug—actsOn→Protein—encodedBy→Gene—hasPhenotype→HPO

## 3. Pipeline
- ETL → RDF/TTL → SHACL → SPARQL (R1)
- Demo: r1_paths output (toy)
- Triple counts & predicate usage

## 4. Evaluation Setup
- Gold: SIDER per-drug AE list (demo-sized)
- Baseline: string-match
- Metrics: P@10 / nDCG@10 / Coverage

## 5. Results
- Toy run: table of top-10 (rule vs baseline) + metrics
- What works vs gaps (HPO↔MedDRA, CPIC β)

## 6. Limitations & Next
- Toy mappings; limited coverage
- Next: real IDs, mapping expansion, CPIC prior weighting, runtime tests

## 7. Reproducibility
- Repo structure, `reports/midterm_2025-10-30` artifacts
- Data source log, licenses, hashes (WIP)
