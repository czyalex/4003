#!/usr/bin/env python3
import os, sys, requests, certifi
from pathlib import Path

URLS = {
    "meddra_all_se.tsv.gz": "https://sideeffects.embl.de/download/meddra_all_se.tsv.gz",
    "drug_names.tsv.gz":    "https://sideeffects.embl.de/download/drug_names.tsv.gz",  # 可选
}

OUTDIR = Path("data/raw/sider")

def download(url: str, out: Path, allow_insecure_fallback: bool = True):
    out.parent.mkdir(parents=True, exist_ok=True)
    def _save(resp):
        resp.raise_for_status()
        with open(out, "wb") as f:
            for chunk in resp.iter_content(1<<20):
                if chunk: f.write(chunk)
        print(f"OK  {out}  ({out.stat().st_size} bytes)")

    try:
        print(f"Downloading (verified): {url}")
        with requests.get(url, stream=True, timeout=60, verify=certifi.where()) as r:
            _save(r)
            return
    except requests.exceptions.SSLError as e:
        if not allow_insecure_fallback:
            raise
        print(f"[SSL error] {e}; trying insecure fallback (verify=False)...")
        with requests.get(url, stream=True, timeout=60, verify=False) as r:
            _save(r)

def main():
    # 如果在代理环境，requests 会自动读取 HTTPS_PROXY/http_proxy 环境变量；也可手动设置：
    # os.environ["HTTPS_PROXY"] = "http://user:pass@proxy.host:port"
    for name, url in URLS.items():
        out = OUTDIR / name
        download(url, out, allow_insecure_fallback=True)

if __name__ == "__main__":
    sys.exit(main())
