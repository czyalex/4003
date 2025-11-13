#!/usr/bin/env python3
"""Toy baseline that 'ranks' AEs by simple string overlap frequency.

Inputs:
- A list of AE strings per drug (e.g., parsed from SIDER).
- A list of HPO labels (or MedDRA PTs mapped to HPO).

This is a stub—replace with your real parser & frequency counts.
"""
from collections import Counter

def rank_ae_by_overlap(ae_strings, vocab):
    counts = Counter()
    vocab_lower = [v.lower() for v in vocab]
    for s in ae_strings:
        low = s.lower()
        for v in vocab_lower:
            if v in low:
                counts[v] += 1
    return [term for term, _ in counts.most_common()]
