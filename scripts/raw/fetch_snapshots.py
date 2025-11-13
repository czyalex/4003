#!/usr/bin/env python3
import sys, csv, hashlib, os, datetime, requests
from pathlib import Path
import yaml

ROOT = Path("."); RAW = ROOT/"data"/"raw"; DOCS = ROOT/"docs"
LOG = DOCS/"data_sources.tsv"

def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""): h.update(chunk)
    return h.hexdigest()

def ensure_log():
    LOG.parent.mkdir(parents=True, exist_ok=True)
    if not LOG.exists():
        LOG.write_text("source\tname\turl\taccessed_on\tversion\tlicense\thash_sha256\tnotes\n", encoding="utf-8")

def append_log(source, name, url, sha, notes="snapshot"):
    accessed = datetime.date.today().isoformat()
    with open(LOG, "a", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([source, name, url, accessed, "", "", sha, notes])

def main(cfg_path):
    ensure_log()
    cfg = yaml.safe_load(open(cfg_path, "r", encoding="utf-8"))
    day = datetime.date.today().isoformat()
    for source, files in cfg.items():
        dest_dir = RAW/source/day; dest_dir.mkdir(parents=True, exist_ok=True)
        print(f"== {source} => {dest_dir}")
        for item in files or []:
            url = item.get("url"); name = item.get("name") or os.path.basename(url or "file")
            if not url or url.startswith("<PUT_"):
                print(f"   - SKIP (fill url later): {name}"); continue
            out = dest_dir/name; print(f"   - GET {url} -> {out}")
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(out, "wb") as f:
                    for chunk in r.iter_content(8192):
                        if chunk: f.write(chunk)
            sha = sha256_of(out); append_log(source, name, url, sha)
            print(f"     SHA256: {sha}")
    print("Done.")

if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv)>1 else "configs/raw_sources.yml"
    main(cfg)
