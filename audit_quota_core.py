"""Independent exact checks of the color-capped core and its completion.

This audit builds rank tables independently of the core implementation and
compares with every placement, not only with the old bounded-core solver.
The stress fixtures deliberately cross the cap and stopping thresholds.
"""
from functools import lru_cache
from itertools import combinations, product
from random import Random

from quota_fault_core import quota_core, solve, solve_control


def mask(s):
    return sum(1 << e for e in s)


def binary_rank_table(columns):
    table = []
    for bits in range(1 << len(columns)):
        basis = {}
        for i, v in enumerate(columns):
            if not (bits >> i & 1):
                continue
            while v:
                p = v.bit_length() - 1
                if p in basis:
                    v ^= basis[p]
                else:
                    basis[p] = v
                    break
        table.append(len(basis))
    return table


def robust_table(ranks, n, k):
    # Repeated worst-case deletions independently check every failure <=k.
    current = ranks[:]
    for _ in range(k):
        current = [min([current[s]] +
                       [current[s ^ (1 << e)] for e in range(n) if s >> e & 1])
                   for s in range(1 << n)]
    return [r == ranks[-1] for r in current]


def check_table(ranks, groups, k, b, shortcut=False):
    n = (len(ranks) - 1).bit_length()
    group_masks = [mask(g) for g in groups]
    robust = robust_table(ranks, n, k)
    true = min((s.bit_count() for s in range(1 << n)
                if robust[s] and all((s & g).bit_count() >= b for g in group_masks)),
               default=None)
    rank = lambda s: ranks[mask(s)]
    placement, core, stats = solve(n, rank, groups, k, b, shortcut)
    assert (None if placement is None else len(placement)) == true, (ranks, groups, k, b, placement, true)
    if placement is not None:
        p = mask(placement)
        assert robust[p] and all((p & g).bit_count() >= b for g in group_masks)
    # Strong size-preserving objective property on the whole ground set.
    if all(len(g) >= b for g in groups):
        q = (k + 1) * ranks[-1]
        best_full, best_core = {}, {}
        cm = mask(core)
        for s in range(1 << n):
            size = s.bit_count()
            if size > q or not robust[s]:
                continue
            value = size + sum(max(0, b - (s & g).bit_count()) for g in group_masks)
            best_full[size] = min(best_full.get(size, value), value)
            if not s & ~cm:
                best_core[size] = min(best_core.get(size, value), value)
        # Loops can be removed, but preserving size is promised for nonloop cores.
        nonloops = sum(1 << i for i in range(n) if ranks[1 << i])
        for size in list(best_full):
            values = [size + sum(max(0, b - (s & g).bit_count()) for g in group_masks)
                      for s in range(1 << n)
                      if s.bit_count() == size and robust[s] and not s & ~nonloops]
            if values:
                assert best_core.get(size) == min(values), (groups, k, b, size, core)
    return stats


def control_table(rows):
    # Independent row-by-row matching DP; achievable right-state subsets.
    n = len(rows)
    reachable = {0}
    for row in rows:
        new = set(reachable)
        for s in reachable:
            for e in range(n):
                if row >> e & 1 and not s >> e & 1:
                    new.add(s | (1 << e))
        reachable = new
    full = (1 << n) - 1
    matchrank = [max((t.bit_count() for t in reachable if not t & ~s), default=0)
                 for s in range(1 << n)]
    dualrank = [s.bit_count() + matchrank[full ^ s] - matchrank[full]
                for s in range(1 << n)]
    # Independent source SCC detection by repeated graph traversal.
    reach = []
    for start in range(n):
        seen = {start}
        todo = [start]
        for v in todo:
            for w in range(n):
                if rows[v] >> w & 1 and w not in seen:
                    seen.add(w)
                    todo.append(w)
        reach.append(seen)
    components = {frozenset(j for j in range(n) if j in reach[i] and i in reach[j])
                  for i in range(n)}
    sources = [c for c in components
               if not any(rows[i] >> j & 1 for i in range(n) if i not in c for j in c)]
    return dualrank, sources


def main():
    totals = dict(nodes=0, threshold=0, color_rejected=0, rank_rejected=0)
    checked = 0
    rng = Random(20261106)
    # Every binary rank-two representation through four ground elements,
    # with all assignments to two disjoint groups or the ungrouped color.
    for n in range(1, 5):
        for columns in product(range(4), repeat=n):
            ranks = binary_rank_table(columns)
            for colors in product(range(3), repeat=n):
                groups = [frozenset(i for i, c in enumerate(colors) if c == j)
                          for j in (1, 2)]
                groups = [g for g in groups if g]
                stats = check_table(ranks, groups, 1, 1)
                for key in totals:
                    totals[key] += stats.get(key, 0)
                checked += 1
    # Higher rank, loops, larger quotas, and multiple failures.
    for _ in range(400):
        n = rng.randrange(1, 11)
        r = rng.randrange(1, 5)
        ranks = binary_rank_table([rng.randrange(1 << r) for _ in range(n)])
        colors = [rng.randrange(4) for _ in range(n)]
        groups = [frozenset(i for i, c in enumerate(colors) if c == j)
                  for j in (1, 2, 3)]
        groups = [g for g in groups if g]
        stats = check_table(ranks, groups, rng.randrange(3), rng.randrange(1, 4))
        for key in totals:
            totals[key] += stats.get(key, 0)
        checked += 1
    # Every directed graph through three vertices: exact actual actuator optimum.
    graph_count = 0
    for n in range(1, 4):
        for rows in product(range(1 << n), repeat=n):
            ranks, sources = control_table(rows)
            for k in (0, 1, 2):
                check_table(ranks, sources, k, k + 1)
                answer, _, _ = solve_control(rows, k)
                feasible = robust_table(ranks, n, k)
                optimum = min((s.bit_count() for s in range(1 << n)
                               if feasible[s] and all(sum(s >> i & 1 for i in c) >= k + 1
                                                      for c in sources)), default=None)
                assert (None if answer is None else len(answer)) == optimum
                graph_count += 1
    for _ in range(80):
        n = rng.randrange(4, 9)
        density = rng.choice((0.2, 0.4, 0.7))
        rows = tuple(sum(1 << j for j in range(n) if rng.random() < density)
                     for _ in range(n))
        ranks, sources = control_table(rows)
        for k in (1, 2):
            check_table(ranks, sources, k, k + 1)
            answer, _, _ = solve_control(rows, k)
            feasible = robust_table(ranks, n, k)
            optimum = min((s.bit_count() for s in range(1 << n)
                           if feasible[s] and all(sum(s >> i & 1 for i in c) >= k + 1
                                                  for c in sources)), default=None)
            assert (None if answer is None else len(answer)) == optimum
            graph_count += 1
    # Arbitrarily many source SCCs need not enlarge the matching parameter.
    rows = [6, 1, 1]
    for _ in range(20):
        a = len(rows)
        rows += [1 << (a + 1), 1 << a]
    placement, core, _ = solve_control(rows)
    assert len(placement) == 42 and len(core) == 2
    print("Exact small matroid/quota comparisons:", checked)
    print("Independent directed-control comparisons:", graph_count)
    print("Fixed-deficiency fixture: 43 states, 21 source SCCs, optimum 42, core 2")
    print("Compression branches:", totals)
    threshold_stress()


def threshold_stress():
    # Rank 1: cap q=2 and total threshold 6, with loops reserved for completion.
    ranks = [int(bool(s & ((1 << 14) - 1))) for s in range(1 << 16)]
    groups = [frozenset(range(i, i + 2)) for i in range(0, 16, 2)]
    stats = check_table(ranks, groups, 1, 2)
    assert stats["threshold"]
    # Rank 2 uniform: c=7, threshold=35.  Cover all cap, rank-rejection and
    # threshold branches via parallel classes and independently enumerate
    # every robust core of size <=4 in both original and compressed ground.
    for seed in range(3):
        n = 44
        directions = list(range(n))
        if seed == 1:
            directions = [i // 2 if i < 12 else i for i in range(n)]
        elif seed == 2:
            directions = [i // 3 if i < 18 else i for i in range(n)]
        groups = [frozenset(range(0, 10))] + [frozenset({i, i + 1}) for i in range(10, 42, 2)]
        # Last two vertices are ungrouped, exercising its special color.
        rank = lambda s: min(2, len({directions[e] for e in s}))
        core, stat = quota_core(n, rank, groups, 1, shortcut=False)
        def optimum(ground):
            best = None
            for z in range(5):
                for c in combinations(sorted(ground), z):
                    if any(rank(set(c) - {e}) != 2 for e in c) or rank(set(c)) != 2:
                        continue
                    value = z + sum(max(0, 2 - len(set(c) & g)) for g in groups)
                    best = value if best is None else min(best, value)
            return best
        assert optimum(range(n)) == optimum(core)
        if seed == 0:
            assert stat["threshold"] and stat["color_rejected"]
        print("Threshold fixture", seed, "core", len(core), "stats", stat)
    print("All quota-core checks passed.")


if __name__ == "__main__":
    main()
