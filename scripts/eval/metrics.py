#!/usr/bin/env python3
import math

def p_at_k(ranked_ids, positive_ids, k=10):
    hits = sum(1 for x in ranked_ids[:k] if x in positive_ids)
    return hits / float(k)

def dcg(ranked_ids, positive_ids, k=10):
    s = 0.0
    for i, x in enumerate(ranked_ids[:k], start=1):
        rel = 1.0 if x in positive_ids else 0.0
        s += (2**rel - 1) / math.log2(i + 1)
    return s

def ndcg(ranked_ids, positive_ids, k=10):
    ideal = dcg(sorted(positive_ids, key=lambda _: True)[:k], positive_ids, k)
    if ideal == 0:
        return 0.0
    return dcg(ranked_ids, positive_ids, k) / ideal
