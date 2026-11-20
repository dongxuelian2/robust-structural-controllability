"""Exact small-instance audit with two mandatory singleton resources.

Importing this module or its helper performs no enumeration.
"""
from itertools import combinations, permutations
from random import Random
from audit_bicircular_pathcover import robust


def anchored_pc(n, edges, s, t):
    """Exact path cover with s,t endpoints, allowing singleton paths."""
    edge_set = {frozenset(e) for e in edges}
    allowed = set()
    for k in range(1, n + 1):
        for path in permutations(range(n), k):
            if s in path[1:-1] or t in path[1:-1]:
                continue
            if all(frozenset((path[j], path[j + 1])) in edge_set
                   for j in range(k - 1)):
                allowed.add(sum(1 << v for v in path))
    best = [n] * (1 << n)
    best[0] = 0
    for mask in range(1, 1 << n):
        for subset in allowed:
            if subset & mask == subset:
                best[mask] = min(best[mask], 1 + best[mask ^ subset])
    return best[-1]


def main():
    tested = 0
    rng = Random(21)
    for n in range(2, 6):
        pairs = list(combinations(range(n), 2))
        masks = (range(1 << len(pairs)) if n <= 4 else
                 [rng.randrange(1 << len(pairs)) for _ in range(30)])
        for mask in masks:
            edges = [e for j, e in enumerate(pairs) if mask >> j & 1]
            if any(not any(v in e for e in edges) for v in range(n)):
                continue
            for s, t in combinations(range(n), 2):
                paths = anchored_pc(n, edges, s, t)
                required = [(s, s), (t, t)]
                optional = edges + [(v, v) for v in range(n) if v not in (s, t)]
                optimum = None
                for k in range(n + 1, n + paths + 1):
                    if any(robust(n, list(selected) + required)
                           for selected in combinations(optional, k - 2)):
                        optimum = k
                        break
                assert optimum == n + paths, (n, edges, s, t, paths, optimum)
                tested += 1
    print(f'PASS: {tested} two-anchor path-cover instances.')


if __name__ == '__main__':
    main()
