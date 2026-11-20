"""Exact path-cover identity audit without singleton completion."""
from itertools import combinations
from random import Random
from audit_bicircular_pathcover import pc, robust


def main():
    tested = 0
    rng = Random(12)
    for n in range(1, 8):
        pairs = list(combinations(range(n), 2))
        masks = (range(1 << len(pairs)) if n <= 5 else
                 [rng.randrange(1 << len(pairs)) for _ in range(200)])
        for mask in masks:
            edges = [e for j, e in enumerate(pairs) if mask >> j & 1]
            if not robust(n, edges):
                continue
            paths = pc(n, edges)
            optimum = None
            for k in range(n + 1, n + paths + 1):
                if any(robust(n, selected) for selected in combinations(edges, k)):
                    optimum = k
                    break
            assert optimum == n + paths, (n, edges, paths, optimum)
            tested += 1
    print(f'PASS: {tested} uncompleted path-cover instances.')


if __name__ == '__main__':
    main()
