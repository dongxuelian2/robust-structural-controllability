"""Independent brute-force check for the truncated-partition mixed-fault DP.

All selected subsets and all admissible failures are compared with the
order-statistic threshold. The optimal values are checked independently.
"""
from itertools import combinations
from random import Random

from truncated_partition_fault import solve


def make_oracle(class_ids, r):
    def independent(elements):
        cls = [class_ids[e] for e in elements]
        return (len(elements) <= r and -1 not in cls
                and len(cls) == len(set(cls)))

    def rank(elements):
        return min(r, len({class_ids[e] for e in elements
                           if class_ids[e] != -1}))
    return independent, rank


def brute(weights, rank, r, k, vulnerable):
    n = len(weights)
    optimum = witness = None
    for size in range(n + 1):
        for combo in combinations(range(n), size):
            chosen = set(combo)
            cost = sum(weights[x] for x in chosen)
            if optimum is not None and cost > optimum:
                continue
            can_fail = sorted(chosen & vulnerable)
            feasible = all(
                rank(chosen.difference(failure)) == r
                for m in range(min(k, len(can_fail)) + 1)
                for failure in combinations(can_fail, m)
            )
            if feasible and (optimum is None or cost < optimum):
                optimum, witness = cost, frozenset(chosen)
    return optimum, witness


def threshold(chosen, class_ids, r, k, vulnerable):
    if r == 0:
        return True
    effective = min(k, len(vulnerable))
    p = max(class_ids, default=-1) + 1
    d = [0] * p
    for x in chosen:
        class_id = class_ids[x]
        if class_id < 0:
            continue
        if x in vulnerable:
            d[class_id] += 1
        else:
            d[class_id] = effective + 1
    d = [min(c, effective + 1) for c in d]
    ordered = sum(sorted(d)[:p-r+1]) >= effective + 1
    dual = max(
        sum(min(t, c) for c in d) - (r-1)*t
        for t in range(effective+2)
    ) >= effective + 1
    assert ordered == dual
    return ordered


rng = Random(20261116)
fixtures = subsets = 0
for case in range(550):
    p = rng.randrange(1, 7)
    sizes = [rng.randrange(1, 4) for _ in range(p)]
    if sum(sizes) > 11:
        continue
    ids = [i for i, size in enumerate(sizes) for _ in range(size)]
    ids += [-1] * rng.randrange(0, 3)
    r = rng.randrange(0, p + 1)
    n = len(ids)
    weights = tuple(rng.randrange(0, 9) for _ in range(n))
    vulnerable = frozenset(i for i in range(n) if rng.randrange(2))
    k = rng.choice([0, 1, 2, 3, n, n+2])
    independent, rank = make_oracle(ids, r)
    actual, selected, certificate = solve(weights, independent, k, vulnerable)
    optimum, witness = brute(weights, rank, r, k, vulnerable)
    assert actual == optimum, (case, ids, r, k, weights, sorted(vulnerable),
                               actual, optimum, selected, witness)
    if selected is not None:
        assert rank(selected) == r or r == 0
        assert sum(weights[x] for x in selected) == actual
    if r > 0:
        for mask in range(1 << n):
            chosen = frozenset(i for i in range(n) if mask >> i & 1)
            can_fail = sorted(chosen & vulnerable)
            exact = all(rank(chosen.difference(failure)) == r
                        for j in range(min(k, len(can_fail)) + 1)
                        for failure in combinations(can_fail, j))
            assert threshold(chosen, ids, r, k, vulnerable) == exact, (
                case, chosen, exact
            )
            subsets += 1
    fixtures += 1

print(f"PASS: {fixtures} fixtures; {subsets} all-subset threshold checks; "
      "all optima match brute force")
