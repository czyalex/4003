# Midterm demo runner (Windows PowerShell) — run from repo root
Continue = "Stop"
python scripts/etl/build_graph.py
python scripts/make_midterm_artifacts.py
Write-Host "✅ Midterm artifacts created. See .\reports\midterm_YYYY-MM-DD"
