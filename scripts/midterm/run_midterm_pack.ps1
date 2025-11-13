param(
  [switch]$WithPTFreq = $false  # 勾上则用“R1+CPIC+PT频次先验”的排名字脚本（若你已有）
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Step($msg){ "`n=== $msg ===" }

# 0) 准备输出目录
$stamp = (Get-Date).ToString("yyyy-MM-dd_HHmm")
$FORMAL_DIR = (Get-ChildItem .\reports\formal_* -Directory | Sort-Object LastWriteTime -Desc | Select-Object -First 1).FullName
$MID = "reports\midterm_$stamp"
New-Item -ItemType Directory -Force -Path $MID | Out-Null
New-Item -ItemType Directory -Force -Path "$MID\figs" | Out-Null

# 1) 快速体检
Step "Smoke check"
& .\scripts\validate\smoke_check.ps1

# 2) 确保 gene→HPO 覆盖结构∪CPIC 基因；确保 HPO→PT 映射存在
Step "Build gene→HPO (include CPIC genes)"
python scripts/interim/make_gene_hpo_from_hpo.py

if (!(Test-Path data\interim\mappings\hpo_meddra_map.tsv)) {
  Step "Build HPO→PT map (hp.obo + SIDER)"
  python scripts/eval/build_hpo2pt_map.py
}

# 3) 重建图谱 & 排名 & 评测
Step "Build graph"
python scripts/etl/build_graph.py

Step ("Rank (" + ($(if($WithPTFreq){"R1+CPIC+PTfreq"} else {"R1+CPIC"}) ) + ")")
if($WithPTFreq){
  python scripts/rankers/rule_r1_with_cpic_and_ptfreq.py
}else{
  python scripts/rankers/rule_r1_with_cpic.py
}

Step "Evaluate (writes latest formal_*/eval_metrics.csv)"
python scripts/eval/run_eval_formal.py

# 4) 生成展图
Step "Make figures"
python scripts/figs/make_midterm_figs.py

# 5) 生成 Top-10 预览（HPO→PT 表头对齐）
Step "Make top10 preview"
python scripts/eval/make_topk_preview.py

# 6) 收集文件到 midterm 目录
$FORMAL_DIR = (Get-ChildItem .\reports\formal_* -Directory | Sort-Object LastWriteTime -Desc | Select-Object -First 1).FullName
Copy-Item $FORMAL_DIR\eval_metrics.csv $MID -Force
Copy-Item $FORMAL_DIR\debug_*.txt      $MID -Force -ErrorAction SilentlyContinue
Copy-Item $FORMAL_DIR\CHEMBL*_ranked_hpo.csv $MID -Force
Copy-Item $FORMAL_DIR\hpo2pt_coverage.txt $MID -Force -ErrorAction SilentlyContinue

Copy-Item figs\*.png                   "$MID\figs" -Force
Copy-Item data\interim\mappings\*.csv  $MID -Force
Copy-Item data\interim\mappings\*.tsv  $MID -Force
if(Test-Path .\triple_counts.txt){ Copy-Item .\triple_counts.txt $MID -Force }
if(Test-Path .\predicate_usage.csv){ Copy-Item .\predicate_usage.csv $MID -Force }

# 7) 生成一个简短 SUMMARY.txt
$met = Import-Csv (Join-Path $FORMAL_DIR 'eval_metrics.csv')
$lines = @()
$lines += "Midterm summary  ($stamp)"
$lines += "Ranking: " + ($(if($WithPTFreq){"R1 + β·CPIC + λ·PTfreq"} else {"R1 + β·CPIC"}))
$lines += ""
$lines += "Metrics:"
$met | ForEach-Object { $lines += (" - {0}: P@10={1}, nDCG@10={2}" -f $_.drug, $_."P@10", $_."nDCG@10") }
$lines += ""
$lines += "Artifacts:"
$lines += " - eval_metrics.csv, debug_*.txt, CHEMBL*_ranked_hpo.csv"
$lines += " - figs/*.png (p10_bar / ndcg_bar / coverage / targets_per_drug)"
$lines += " - hpo2pt_coverage.txt (if present)"
$lines += " - mappings/*.csv|*.tsv"
$lines | Set-Content -Encoding UTF8 (Join-Path $MID 'SUMMARY.txt')

# 8) 打包 ZIP
$zip = "$MID.zip"
if(Test-Path $zip){ Remove-Item $zip -Force }
Compress-Archive -Path "$MID\*" -DestinationPath $zip -Force
"Done. Midterm pack: $zip"
