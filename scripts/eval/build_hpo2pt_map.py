#!/usr/bin/env python3
"""
Build HPO -> MedDRA PT mapping for AE-KG.

Inputs (matching your 4003 repo):
- gene_hpo.csv: gene_id, hpo_id (list of HPO terms actually used)
- hp.obo: full HPO ontology in OBO format (for labels and synonyms)
- pt_synonyms.tsv: mapping from LLT to PT
    columns: llt, pt
    (we treat 'pt' as both PT code and PT name for now)
- hpo_meddra_overrides.tsv: manual overrides (hpo_id, pt_code)

Outputs:
- hpo_meddra_map.tsv: final mapping per HPO
- coverage txt file with summary stats (JSON text)
"""

import argparse
import json
import re
from collections import defaultdict
from typing import Dict, List, Set

import pandas as pd


# ---------- Text utilities ----------

def normalize_text(s: str) -> str:
    """Normalize text for matching (lowercase, trim, remove simple punctuation, unify arrows)."""
    if not isinstance(s, str):
        return ""
    s = s.lower().strip()
    for ch in [",", ";", ":", "(", ")", "[", "]"]:
        s = s.replace(ch, " ")
    s = s.replace("↑", " increased ")
    s = s.replace("↓", " decreased ")
    s = " ".join(s.split())
    return s


def split_synonyms(s: str) -> List[str]:
    """Split synonyms if needed. Here it mostly returns singletons, but keep it generic."""
    if not isinstance(s, str) or not s.strip():
        return []
    for sep in ["|", ";;"]:
        if sep in s:
            return [p.strip() for p in s.split(sep) if p.strip()]
    return [s.strip()]


# ---------- HPO side: parse hp.obo ----------

def load_hpo_terms_from_obo(obo_path: str) -> Dict[str, Dict[str, List[str]]]:
    """
    Parse hp.obo and return:
        dict[hpo_id] = {"label": name, "synonyms": [syn1, syn2, ...]}
    Only uses [Term] blocks with 'id:' and 'name:'.
    """
    hpo_terms: Dict[str, Dict[str, List[str]]] = {}

    current_id = None
    current_name = None
    current_syns: List[str] = []
    in_term = False

    id_prefix = "id: "
    name_prefix = "name: "
    syn_prefix = "synonym: "

    syn_re = re.compile(r'^synonym: "(.+?)"')

    with open(obo_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")

            if line == "[Term]":
                # flush previous term
                if in_term and current_id and current_name:
                    hpo_terms[current_id] = {
                        "label": current_name,
                        "synonyms": current_syns,
                    }
                # start new term
                in_term = True
                current_id = None
                current_name = None
                current_syns = []
                continue

            if not line or line.startswith("[Typedef]"):
                # end of a [Term] block
                if in_term and current_id and current_name:
                    hpo_terms[current_id] = {
                        "label": current_name,
                        "synonyms": current_syns,
                    }
                in_term = False
                current_id = None
                current_name = None
                current_syns = []
                continue

            if not in_term:
                continue

            if line.startswith(id_prefix):
                current_id = line[len(id_prefix):].strip()
            elif line.startswith(name_prefix):
                current_name = line[len(name_prefix):].strip()
            elif line.startswith(syn_prefix):
                m = syn_re.match(line)
                if m:
                    syn_txt = m.group(1).strip()
                    if syn_txt:
                        current_syns.append(syn_txt)

    # flush last term if file did not end with blank/[Typedef]
    if in_term and current_id and current_name:
        hpo_terms[current_id] = {
            "label": current_name,
            "synonyms": current_syns,
        }

    return hpo_terms


def build_hpo_strings(
    hpo_ids: List[str],
    hpo_terms: Dict[str, Dict[str, List[str]]],
) -> Dict[str, List[str]]:
    """
    For each HPO ID in hpo_ids, collect strings: label + synonyms (if available).
    Returns: dict[hpo_id] -> list of strings.
    """
    hpo_strings: Dict[str, List[str]] = {}

    for hid in hpo_ids:
        info = hpo_terms.get(hid)
        if not info:
            # No entry in hp.obo; keep empty list
            hpo_strings[hid] = []
            continue

        label = info.get("label") or ""
        syns = info.get("synonyms") or []

        strings: List[str] = []
        if label:
            strings.append(label)
        strings.extend(syns)

        # de-duplicate
        seen: Set[str] = set()
        uniq: List[str] = []
        for s in strings:
            s = s.strip()
            if s and s not in seen:
                seen.add(s)
                uniq.append(s)

        hpo_strings[hid] = uniq

    return hpo_strings


# ---------- PT side: build index from pt_synonyms.tsv ----------

def build_pt_index(
    df_pt: pd.DataFrame,
    pt_code_col: str,
    pt_name_col: str,
    llt_col: str,
) -> tuple[Dict[str, Set[str]], Dict[str, str]]:
    """
    Build:
    - norm_to_pt: normalized text (PT name or LLT) -> set of PT codes
    - pt_code_to_name: PT code -> PT name

    In your current file:
        pt_code_col = "pt"
        pt_name_col = "pt"
        llt_col = "llt"
    """
    norm_to_pt: Dict[str, Set[str]] = defaultdict(set)
    pt_code_to_name: Dict[str, str] = {}

    for _, row in df_pt.iterrows():
        code = row.get(pt_code_col)
        name = row.get(pt_name_col)
        llt = row.get(llt_col)

        if not pd.isna(code):
            code = str(code)
        else:
            continue

        # Use PT name as canonical name (may be same as code)
        if not pd.isna(name):
            name = str(name)
        else:
            name = code

        pt_code_to_name[code] = name

        # Index PT name itself
        norm_name = normalize_text(name)
        if norm_name:
            norm_to_pt[norm_name].add(code)

        # Index LLT as synonym
        if not pd.isna(llt):
            llt_str = str(llt)
            for syn in split_synonyms(llt_str):
                norm_llt = normalize_text(syn)
                if norm_llt:
                    norm_to_pt[norm_llt].add(code)

    return norm_to_pt, pt_code_to_name


# ---------- Mapping HPO -> PT candidates ----------

def map_hpo_to_pts(
    hpo_strings: Dict[str, List[str]],
    norm_to_pt: Dict[str, Set[str]],
) -> Dict[str, Set[str]]:
    """Use normalized HPO strings to look up PT candidates."""
    hpo_to_pts: Dict[str, Set[str]] = defaultdict(set)

    for hid, strings in hpo_strings.items():
        for s in strings:
            norm = normalize_text(s)
            if not norm:
                continue
            pt_set = norm_to_pt.get(norm)
            if pt_set:
                hpo_to_pts[hid].update(pt_set)

    return hpo_to_pts


# ---------- Overrides ----------

def load_overrides(
    path: str,
    hpo_id_col: str,
    pt_code_col: str,
) -> Dict[str, str]:
    """Load manual overrides: HPO ID -> PT code."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    overrides: Dict[str, str] = {}
    for _, row in df.iterrows():
        hid = row.get(hpo_id_col)
        code = row.get(pt_code_col)
        if pd.isna(hid) or pd.isna(code):
            continue
        overrides[str(hid)] = str(code)
    return overrides


# ---------- Main ----------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build HPO -> MedDRA PT mapping for AE-KG."
    )

    # HPO / gene_hpo inputs
    parser.add_argument(
        "--gene-hpo",
        required=True,
        help="Path to gene_hpo.csv (with hpo_id column).",
    )
    parser.add_argument(
        "--hpo-id-col",
        default="hpo_id",
        help="Column name for HPO ID in gene_hpo.csv.",
    )
    parser.add_argument(
        "--hpo-obo",
        required=True,
        help="Path to hp.obo (HPO ontology in OBO format).",
    )

    # PT side (LLT -> PT)
    parser.add_argument(
        "--pt-synonyms",
        required=True,
        help="Path to pt_synonyms.tsv (llt, pt).",
    )
    parser.add_argument(
        "--pt-code-col",
        default="pt",
        help="Column name for PT code (in your file: 'pt').",
    )
    parser.add_argument(
        "--pt-name-col",
        default="pt",
        help="Column name for PT name (in your file: also 'pt').",
    )
    parser.add_argument(
        "--llt-col",
        default="llt",
        help="Column name for LLT text (in your file: 'llt').",
    )

    # Overrides
    parser.add_argument(
        "--overrides",
        default=None,
        help="Optional path to hpo_meddra_overrides.tsv.",
    )
    parser.add_argument(
        "--overrides-hpo-col",
        default="hpo_id",
        help="Column name for HPO ID in overrides file.",
    )
    parser.add_argument(
        "--overrides-pt-col",
        default="pt_code",
        help="Column name for PT code in overrides file.",
    )

    # Outputs
    parser.add_argument(
        "--out-tsv",
        required=True,
        help="Output TSV path for final hpo_meddra_map.tsv.",
    )
    parser.add_argument(
        "--coverage-txt",
        default=None,
        help="Optional path to write coverage stats as JSON text.",
    )

    args = parser.parse_args()

    # ---- Load gene_hpo and HPO terms ----
    df_hpo = pd.read_csv(args.gene_hpo, dtype=str)
    if args.hpo_id_col not in df_hpo.columns:
        raise ValueError(
            f"HPO ID column '{args.hpo_id_col}' not found in {args.gene_hpo}. "
            f"Available columns: {df_hpo.columns.tolist()}"
        )

    hpo_ids: List[str] = (
        df_hpo[args.hpo_id_col]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    hpo_terms = load_hpo_terms_from_obo(args.hpo_obo)
    hpo_strings = build_hpo_strings(hpo_ids, hpo_terms)

    # ---- Load PT synonyms and build index ----
    df_pt = pd.read_csv(
        args.pt_synonyms,
        sep=None,
        engine="python",
        dtype=str,
    )
    # Strip BOM and whitespace from column names
    df_pt.columns = [c.strip().lstrip("\ufeff") for c in df_pt.columns]

    if args.pt_code_col not in df_pt.columns or args.pt_name_col not in df_pt.columns:
        raise ValueError(
            f"PT columns '{args.pt_code_col}'/'{args.pt_name_col}' not found in pt_synonyms. "
            f"Available columns: {df_pt.columns.tolist()}"
        )
    if args.llt_col not in df_pt.columns:
        raise ValueError(
            f"LLT column '{args.llt_col}' not found in pt_synonyms. "
            f"Available columns: {df_pt.columns.tolist()}"
        )

    norm_to_pt, pt_code_to_name = build_pt_index(
        df_pt,
        pt_code_col=args.pt_code_col,
        pt_name_col=args.pt_name_col,
        llt_col=args.llt_col,
    )

    # ---- Map HPO -> PT candidates ----
    hpo_to_pts = map_hpo_to_pts(hpo_strings, norm_to_pt)

    # ---- Overrides ----
    overrides: Dict[str, str] = {}
    if args.overrides:
        overrides = load_overrides(
            path=args.overrides,
            hpo_id_col=args.overrides_hpo_col,
            pt_code_col=args.overrides_pt_col,  # type: ignore[attr-defined]
        )

    # ---- Build final table ----
    rows: List[dict] = []
    n_total = len(hpo_strings)
    n_override = n_unique = n_ambiguous = n_no = 0

    for hid, strings in hpo_strings.items():
        candidates = sorted(hpo_to_pts.get(hid, set()))
        cand_names = [pt_code_to_name.get(c, "") for c in candidates]

        method = ""
        chosen_code = ""
        chosen_name = ""

        if hid in overrides:
            chosen_code = overrides[hid]
            chosen_name = pt_code_to_name.get(chosen_code, chosen_code)
            method = "override"
            n_override += 1
        else:
            if not candidates:
                method = "no_match"
                n_no += 1
            elif len(candidates) == 1:
                chosen_code = candidates[0]
                chosen_name = pt_code_to_name.get(chosen_code, chosen_code)
                method = "unique_auto_match"
                n_unique += 1
            else:
                method = "ambiguous"
                n_ambiguous += 1

        rows.append(
            {
                "hpo_id": hid,
                "hpo_strings": "|".join(strings),
                "pt_code": chosen_code,
                "pt_name": chosen_name,
                "method": method,
                "n_candidates": len(candidates),
                "candidate_pt_codes": "|".join(candidates),
                "candidate_pt_names": "|".join(cand_names),
            }
        )

    df_out = pd.DataFrame(rows)
    df_out.to_csv(args.out_tsv, sep="\t", index=False)

    # ---- Coverage stats ----
    stats = {
        "n_total_hpo": n_total,
        "n_override": n_override,
        "n_unique_auto": n_unique,
        "n_ambiguous": n_ambiguous,
        "n_no_match": n_no,
        "covered_hpo": n_total - n_no,
        "coverage_ratio": (n_total - n_no) / n_total if n_total else 0.0,
    }

    txt = json.dumps(stats, indent=2, sort_keys=True)
    print(txt)

    if args.coverage_txt:
        with open(args.coverage_txt, "w", encoding="utf-8") as f:
            f.write(txt + "\n")


if __name__ == "__main__":
    main()
