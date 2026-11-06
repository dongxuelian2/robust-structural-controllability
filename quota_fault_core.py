"""Unit-cost fault-tolerant spanning with disjoint lower quotas.

All selected elements may fail.  The parameter is matroid rank plus the
failure budget, not the number of groups.  Integer elements are 0,...,n-1;
groups are disjoint and each has the common positive quota b (default k+1).
This reference solver favors a short auditable implementation, not speed.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb


def _inputs(n, groups, k, b):
    if type(n) is not int or n < 0 or type(k) is not int or k < 0:
        raise ValueError("n and k must be nonnegative integers")
    if type(b) is not int or b < 1:
        raise ValueError("b must be a positive integer")
    groups = tuple(frozenset(g) for g in groups)
    colors = [0] * n
    for i, group in enumerate(groups, 1):
        for e in group:
            if type(e) is not int or not 0 <= e < n or colors[e]:
                raise ValueError("groups must be disjoint subsets of 0,...,n-1")
            colors[e] = i
    return groups, tuple(colors)


def quota_core(n, rank_oracle, groups, k, b=None, shortcut=True):
    """Return a small robust-core candidate set, with execution statistics.

    The returned set is NOT itself necessarily a feasible quota placement.
    Quota completions may use original elements outside the returned core.
    shortcut=False exercises compression even on small test instances.
    """
    b = k + 1 if b is None else b
    groups, colors = _inputs(n, groups, k, b)

    @lru_cache(None)
    def rank(a):
        return rank_oracle(a)

    ground = frozenset(range(n))
    d = rank(ground)
    q = (k + 1) * d
    stats = dict(nodes=0, threshold=0, color_rejected=0, rank_rejected=0)
    if d == 0:
        return frozenset(), stats
    H = comb(d + k - 1, k)

    @lru_cache(None)
    def visit(a):
        stats["nodes"] += 1
        h = rank(a)
        if h == 0:
            return frozenset()
        cap = (h - 1) * H + q
        threshold = (q + 1) * cap
        if shortcut and len(a) <= threshold:
            return a
        z = []
        counts = Counter()
        for x in sorted(a):
            if counts[colors[x]] == cap:
                stats["color_rejected"] += 1
                continue
            if len(z) < h - 1:
                accept = rank(frozenset(z + [x])) == len(z) + 1
            else:
                accept = all(rank(frozenset(t + (x,))) == h
                             for t in combinations(z, h - 1))
            if not accept:
                stats["rank_rejected"] += 1
                continue
            z.append(x)
            counts[colors[x]] += 1
            if len(z) == threshold:
                stats["threshold"] += 1
                break
        assert len(z) >= h
        result = set(z)
        if h >= 2:
            for t in combinations(z, h - 1):
                tf = frozenset(t)
                child = frozenset(x for x in a if rank(tf | {x}) == h - 1)
                assert rank(child) == h - 1
                result.update(visit(child))
        return frozenset(result)

    nonloops = frozenset(e for e in ground if rank(frozenset({e})))
    return visit(nonloops), stats


def solve(n, rank_oracle, groups, k=1, b=None, shortcut=True):
    """Return (optimum placement or None, candidate core, statistics)."""
    b = k + 1 if b is None else b
    groups, _ = _inputs(n, groups, k, b)
    if any(len(g) < b for g in groups):
        return None, frozenset(), {}
    core, stats = quota_core(n, rank_oracle, groups, k, b, shortcut)
    d = rank_oracle(frozenset(range(n)))
    q = (k + 1) * d
    best = None
    for size in range(min(q, len(core)) + 1):
        for chosen in combinations(sorted(core), size):
            c = frozenset(chosen)
            cost = size + sum(max(0, b - len(c & g)) for g in groups)
            if best is not None and cost >= len(best):
                continue
            if not all(rank_oracle(c - frozenset(f)) == d
                       for f in combinations(chosen, min(k, size))):
                continue
            placement = set(c)
            for g in groups:
                need = max(0, b - len(c & g))
                placement.update(sorted(g - c)[:need])
            assert len(placement) == cost
            best = frozenset(placement)
    return best, core, stats


def solve_control(rows, k=1):
    """Exact general-digraph dedicated-actuator placement under k failures."""
    from bounded_core_solver import matching_number, source_components
    rows = tuple(rows)
    n = len(rows)
    if not n or any(type(row) is not int or row < 0 or row >= 1 << n
                    for row in rows):
        raise ValueError("expected nonempty nonnegative n-bit adjacency rows")
    full = (1 << n) - 1
    nu = matching_number(rows, full)
    groups = [frozenset(i for i in range(n) if c >> i & 1)
              for c in source_components(rows)]

    @lru_cache(None)
    def rank(s):
        mask = sum(1 << e for e in s)
        return len(s) + matching_number(rows, full ^ mask) - nu

    return solve(n, rank, groups, k)
