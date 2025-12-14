## MVP next steps

1. Fill `docs/drugs/` dossiers for simvastatin & carbamazepine (IDs, CPIC link, AE seed terms).
2. Populate `docs/data_sources.tsv` as you download raw data into `data/raw/`.
3. Create/activate env:
   ```bash
   conda env create -f env/environment.yml
   conda activate ae-kg
   pip install -r env/requirements.txt
   ```
4. Make placeholder raw folders:
   ```bash
   python scripts/download_raw.py --drug simvastatin --drug carbamazepine
   ```
5. Build a toy graph from interim mappings:
   ```bash
   python scripts/etl/build_graph.py
   ```
6. Run the R1 SPARQL query against your triple store (or load `rdf/ae_kg.ttl` into a viewer).
