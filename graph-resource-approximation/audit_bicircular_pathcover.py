"""Exact small-instance audit of the unanchored path-cover identity.

Importing this module performs no enumeration. Python 3.10+, standard library.
"""
from itertools import combinations
from random import Random


def pc(n, edges):
    """Minimum spanning path-cover size, using exact subset dynamic programming."""
    size = 1 << n
    ends = [0] * size
    adjacency = [0] * n
    for u, v in edges:
        adjacency[u] |= 1 << v
        adjacency[v] |= 1 << u
    for v in range(n):
        ends[1 << v] = 1 << v
    for mask in range(1, size):
        for v in range(n):
            if ends[mask] >> v & 1:
                available = adjacency[v] & ~mask
                while available:
                    bit = available & -available
                    available -= bit
                    ends[mask | bit] |= bit
    best = [n] * size
    best[0] = 0
    for mask in range(1, size):
        subset = mask
        while subset:
            if ends[subset]:
                best[mask] = min(best[mask], 1 + best[mask ^ subset])
            subset = (subset - 1) & mask
    return best[-1]


def robust(n, edges):
    """Graph criterion for robust incidence assignment; a loop is (v,v)."""
    degree = [0] * n
    parent = list(range(n))

    def root(v):
        while parent[v] != v:
            v = parent[v]
        return v

    for u, v in edges:
        degree[u] += 1
        degree[v] += 1
        parent[root(u)] = root(v)
    if min(degree) < 2:
        return False
    vertices = {}
    degree_sum = {}
    for v in range(n):
        r = root(v)
        vertices[r] = vertices.get(r, 0) + 1
        degree_sum[r] = degree_sum.get(r, 0) + degree[v]
    return all(degree_sum[r] > 2 * vertices[r] for r in vertices)


def main():
    tested = 0
    rng = Random(12)
    for n in range(2, 7):
        pairs = list(combinations(range(n), 2))
        masks = (range(1 << len(pairs)) if n <= 5 else
                 [rng.randrange(1 << len(pairs)) for _ in range(100)])
        for mask in masks:
            edges = [e for j, e in enumerate(pairs) if mask >> j & 1]
            if any(not any(v in e for e in edges) for v in range(n)):
                continue
            paths = pc(n, edges)
            resources = edges + [(v, v) for v in range(n)]
            optimum = None
            for k in range(n + 1, n + paths + 1):
                if any(robust(n, selected) for selected in combinations(resources, k)):
                    optimum = k
                    break
            assert optimum == n + paths, (n, edges, paths, optimum)
            tested += 1
    print(f'PASS: {tested} singleton-completed path-cover instances.')


if __name__ == '__main__':
    main()
