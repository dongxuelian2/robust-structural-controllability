"""Exhaustively audit the pendant reduction and its 8/7 ratio calculation.

This tests the reduction on every eligible returned path cover, not an
implementation of the cited Berman-Karpinski approximation algorithm.
"""
from itertools import combinations


def linear_forest(n, edges):
    degree = [0] * n
    parent = list(range(n))

    def root(v):
        while parent[v] != v:
            v = parent[v]
        return v

    for u, v in edges:
        degree[u] += 1
        degree[v] += 1
        if degree[u] > 2 or degree[v] > 2 or root(u) == root(v):
            return False
        parent[root(u)] = root(v)
    return True


def normalize(n, edges, s, t):
    chosen = set(edges)
    for v, leaf in [(s, n), (t, n + 1)]:
        edge = (v, leaf)
        if edge in chosen:
            continue
        incident = [e for e in chosen if v in e]
        if len(incident) == 2:
            chosen.remove(min(incident))
        chosen.add(edge)
    return chosen


def main():
    instances = forests = eligible = 0
    for n in range(2, 5):
        base = list(combinations(range(n), 2))
        for mask in range(1 << len(base)):
            graph = [e for i, e in enumerate(base) if mask >> i & 1]
            if any(not any(v in e for e in graph) for v in range(n)):
                continue
            for s, t in combinations(range(n), 2):
                augmented = graph + [(s, n), (t, n + 1)]
                covers = []
                anchored_best = 0
                for subset in range(1 << len(graph)):
                    chosen = [e for i, e in enumerate(graph) if subset >> i & 1]
                    if (linear_forest(n, chosen) and
                            all(sum(v in e for e in chosen) <= 1 for v in (s, t))):
                        anchored_best = max(anchored_best, len(chosen))
                for subset in range(1 << len(augmented)):
                    chosen = [e for i, e in enumerate(augmented) if subset >> i & 1]
                    if linear_forest(n + 2, chosen):
                        covers.append(chosen)
                optimum = max(map(len, covers))
                assert optimum == anchored_best + 2
                p = n - anchored_best
                for chosen in covers:
                    normalized = normalize(n, chosen, s, t)
                    assert linear_forest(n + 2, normalized)
                    assert len(normalized) >= len(chosen)
                    recovered = {e for e in normalized if e[1] < n}
                    assert len(recovered) == len(normalized) - 2
                    assert linear_forest(n, recovered)
                    assert all(sum(v in e for e in recovered) <= 1 for v in (s, t))
                    q = n - len(recovered)
                    assert q == n + 2 - len(normalized)
                    if 7 * len(chosen) >= 6 * optimum:
                        assert 7 * (n + q) <= 8 * (n + p)
                        eligible += 1
                    forests += 1
                instances += 1
    print(f'PASS: {instances} pendant-anchor instances; {forests} path covers; '
          f'{eligible} covers satisfying the cited 6/7 edge guarantee.')


if __name__ == '__main__':
    main()
