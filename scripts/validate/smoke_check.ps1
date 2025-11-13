# scripts/validate/smoke_check.ps1  —— 修正 ${var}: 用法 + 更稳的检查
param()

function Test-Regex($s, $pat){ return [bool]([regex]$pat).IsMatch($s) }

$fail = @()

# 1) drug_targets.csv
$dtPath = "data/interim/mappings/drug_targets.csv"
if(!(Test-Path $dtPath)){ $fail += "MISSING: ${dtPath}" }
else{
  $dt = Import-Csv $dtPath
  if($dt.Count -eq 0){ $fail += "EMPTY: ${dtPath}" }
  else{
    $uniPat = '(^[A-NR-Z][0-9]{5}$)|(^[OPQ][0-9][A-Z0-9]{3}[0-9]$)'
    $bad = $dt | Where-Object { -not (Test-Regex $_.protein_id $uniPat) }
    if($bad){
      $badShow = ($bad | Select-Object -First 3 | ForEach-Object { $_.protein_id } | Sort-Object -Unique) -join ', '
      $fail += "BAD_UNIPROT in ${dtPath}: ${badShow}"
    }
    $simvaOK = $dt | Where-Object { $_.drug_id -eq 'CHEMBL1064' -and $_.protein_id -eq 'P04035' }
    if(-not $simvaOK){ $fail += "SIMVASTATIN missing P04035(HMGCR) in ${dtPath}" }
  }
}

# 2) protein_gene.csv covers all proteins in drug_targets.csv
$pgPath = "data/interim/mappings/protein_gene.csv"
if(!(Test-Path $pgPath)){ $fail += "MISSING: ${pgPath}" }
else{
  $pg = Import-Csv $pgPath
  if($pg.Count -eq 0){ $fail += "EMPTY: ${pgPath}" }
  else{
    $dtProts = (Import-Csv $dtPath | Select-Object -Expand protein_id | Where-Object { $_ } | Sort-Object -Unique)
    $pgProts = ($pg | Select-Object -Expand protein_id | Where-Object { $_ } | Sort-Object -Unique)
    $missing = @()
    foreach($p in $dtProts){ if($p -notin $pgProts){ $missing += $p } }
    if($missing){ $fail += "PROTEIN->GENE mapping missing for: " + (($missing | Select-Object -First 10) -join ', ') }
  }
}

# 3) gene_hpo.csv — HPO format
$ghPath = "data/interim/mappings/gene_hpo.csv"
if(!(Test-Path $ghPath)){ $fail += "MISSING: ${ghPath}" }
else{
  $gh = Import-Csv $ghPath
  $badHpo = $gh | Where-Object { -not (Test-Regex $_.hpo_id '^HP:\d{7}$') }
  if($badHpo){
    $badShow = ($badHpo | Select-Object -First 3 | ForEach-Object {$_.hpo_id}) -join ', '
    $fail += "BAD_HPO in ${ghPath}: ${badShow}"
  }
}

# 4) hpo_meddra_map.tsv — ensure HP:0003198 -> Myopathy
$mapPath = "data/interim/mappings/hpo_meddra_map.tsv"
if(!(Test-Path $mapPath)){ $fail += "MISSING: ${mapPath}" }
else{
  $map = Import-Csv $mapPath -Delimiter "`t"
  $hit = $map | Where-Object { $_.hpo_id -eq 'HP:0003198' -and $_.meddra_pt -eq 'Myopathy' }
  if(-not $hit){ $fail += "HPO→PT map missing HP:0003198 -> Myopathy (for simvastatin scoring)" }
}

# 5) drug_gene_pgx.csv — columns & numeric weights
$pgxPath = "data/interim/mappings/drug_gene_pgx.csv"
if(!(Test-Path $pgxPath)){ $fail += "MISSING: ${pgxPath}" }
else{
  $pgx = Import-Csv $pgxPath
  $need = @('drug_id','gene_id','cpic_weight')
  foreach($c in $need){
    if(-not ($pgx[0].PSObject.Properties.Name -contains $c)){ $fail += "MISSING_COL ${c} in ${pgxPath}" }
  }
  $badW = $pgx | Where-Object { -not [double]::TryParse(($_.cpic_weight -replace ',','.'), [ref]([double]0)) }
  if($badW){
    $wShow = ($badW | Select-Object -First 3 | ForEach-Object { $_.cpic_weight }) -join ', '
    $fail += "BAD_WEIGHT in ${pgxPath}: ${wShow}"
  }
}

if($fail.Count){
  "❌ CHECKS FAILED:"; $fail | ForEach-Object { " - $_" }
  exit 1
}else{
  "✅ ALL CHECKS PASSED"
  "— drug_targets.csv —"; (Import-Csv $dtPath | Select-Object -First 5) | Format-Table -AutoSize
  "— protein_gene.csv —"; (Import-Csv $pgPath | Select-Object -First 5) | Format-Table -AutoSize
  "— hpo_meddra_map.tsv —"; (Import-Csv $mapPath -Delimiter "`t" | Select-Object -First 5) | Format-Table -AutoSize
}
