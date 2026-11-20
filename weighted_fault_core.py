"""Exact FPT weighted fault-tolerant spanning via a bounded weighted core.

The input rank oracle takes a frozenset of integer element IDs. Costs are
nonnegative integers. Only elements in `vulnerable` can fail (default: all).
This deliberately small reference implementation is not intended for large
parameters: its bound is a function of rank and the failure budget.
"""
from functools import lru_cache
from itertools import combinations
from math import comb


def weighted_core(weights, rank_oracle, k, vulnerable=None):
    weights = tuple(weights)
    if type(k) is not int or k < 0:
        raise ValueError("k must be a nonnegative integer")
    if any(type(w) is not int or w < 0 for w in weights):
        raise ValueError("weights must be nonnegative integers")
    ground = frozenset(range(len(weights)))
    vulnerable = ground if vulnerable is None else frozenset(vulnerable)
    if not vulnerable <= ground:
        raise ValueError("vulnerable elements must belong to the ground set")

    @lru_cache(None)
    def rank(a):
        return rank_oracle(a)

    r = rank(ground)
    q = (k + 1) * r
    stats = dict(nodes=0, capped=0, rejected=0, parallel_pruned=0)
    if r == 0:
        return frozenset(), stats
    # A deficient hyperplane has an (r-1)-element basis disjoint from
    # at most k deleted points. Distinct pairs cross-intersect, so the
    # Bollobas set-pairs bound applies (see the accompanying proof).
    hyperplanes = comb(r + k - 1, k)

    @lru_cache(None)
    def visit(a):
        stats["nodes"] += 1
        h = rank(a)
        if h == 0:
            return frozenset()
        ordered = sorted(a, key=lambda e: (weights[e], e))
        if h == 1:
            stats["parallel_pruned"] += max(0, len(a) - q)
            return frozenset(ordered[:q])
        threshold = (h - 1) * hyperplanes + q
        if len(a) <= threshold:
            return a
        z = []
        for x in ordered:
            if len(z) < h - 1:
                accept = rank(frozenset(z + [x])) == len(z) + 1
            else:
                accept = all(rank(frozenset(s + (x,))) == h
                             for s in combinations(z, h - 1))
            if accept:
                z.append(x)
                if len(z) == threshold:
                    stats["capped"] += 1
                    break
            else:
                stats["rejected"] += 1
        assert len(z) >= h
        result = set(z)
        for s in combinations(z, h - 1):
            sf = frozenset(s)
            child = frozenset(x for x in a if rank(sf | {x}) == h - 1)
            assert rank(child) == h - 1
            result.update(visit(child))
        return frozenset(result)

    nonloops = frozenset(x for x in ground if rank(frozenset({x})) > 0)
    core = visit(nonloops & vulnerable) | visit(nonloops - vulnerable)
    return core, stats


def solve(weights, rank_oracle, k, vulnerable=None):
    """Return (minimum cost, selected frozenset, core, stats), or None cost."""
    weights = tuple(weights)
    if vulnerable is not None:
        vulnerable = frozenset(vulnerable)
    core, stats = weighted_core(weights, rank_oracle, k, vulnerable)
    ground = frozenset(range(len(weights)))
    vulnerable = ground if vulnerable is None else frozenset(vulnerable)
    r = rank_oracle(ground)
    q = (k + 1) * r
    if r == 1:
        # A rank-one set survives iff it contains a safe nonloop or at
        # least k+1 vulnerable nonloops. The core already discards loops
        # and retains the cheapest q elements of each status class.
        safe = sorted(core - vulnerable, key=lambda e: (weights[e], e))
        exposed = sorted(core & vulnerable, key=lambda e: (weights[e], e))
        candidates = [frozenset(safe[:1])] if safe else []
        if len(exposed) >= q:
            candidates.append(frozenset(exposed[:q]))
        if not candidates:
            return None, None, core, stats
        witness = min(candidates, key=lambda s: (
            sum(weights[e] for e in s), len(s), tuple(sorted(s))))
        return sum(weights[e] for e in witness), witness, core, stats
    best = None
    witness = None
    for size in range(min(q, len(core)) + 1):
        for candidate in combinations(sorted(core), size):
            cost = sum(weights[e] for e in candidate)
            if best is not None and cost >= best:
                continue
            b = frozenset(candidate)
            v = sorted(b & vulnerable)
            # Checking maximal-size allowed failures suffices by monotonicity.
            if all(rank_oracle(b - frozenset(f)) == r
                   for f in combinations(v, min(k, len(v)))):
                best, witness = cost, b
    return best, witness, core, stats

