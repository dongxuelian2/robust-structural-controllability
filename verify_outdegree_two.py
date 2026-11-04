"""Exact subset-matching audit of the outdegree-two theorem; no dependencies."""
from itertools import product
from random import Random
from outdegree_two_solver import solve, tree_blocks


def matching_subsets(rows):
    reachable = {0}
    for row in rows:
        previous = tuple(reachable)
        for mask in previous:
            available = row & ~mask
            while available:
                bit = available & -available
                available ^= bit
                reachable.add(mask | bit)
    return reachable


def independent_sources(rows):
    n = len(rows)
    reach = []
    for start in range(n):
        visited, queue = {start}, [start]
        for v in queue:
            for w in range(n):
                if rows[v] >> w & 1 and w not in visited:
                    visited.add(w)
                    queue.append(w)
        reach.append(visited)
    components, unseen = [], set(range(n))
    while unseen:
        v = min(unseen)
        component = {w for w in unseen if w in reach[v] and v in reach[w]}
        unseen.difference_update(component)
        if not any(rows[w] >> z & 1 for w in range(n) if w not in component
                   for z in component):
            components.append(sum(1 << w for w in component))
    return components


def brute_robust(rows):
    n = len(rows)
    full = (1 << n) - 1
    matchable = matching_subsets(rows)
    sources = independent_sources(rows)
    controllable = [full ^ s in matchable and all(s & c for c in sources)
                    for s in range(full + 1)]
    robust = []
    for s in range(full + 1):
        if controllable[s] and all(controllable[s ^ (1 << v)]
                                   for v in range(n) if s >> v & 1):
            robust.append(s)
    return matchable, robust


def check(rows):
    n = len(rows)
    full = (1 << n) - 1
    matchable, robust = brute_robust(rows)
    blocks = tree_blocks(rows)
    for subset in range(full + 1):
        predicted = all((subset & c).bit_count() < c.bit_count() for c in blocks)
        assert (subset in matchable) == predicted, (rows, subset, 'decomposition')
    result = solve(iter(rows))
    optimum = min((s.bit_count() for s in robust), default=None)
    assert (None if result is None else len(result)) == optimum, (rows, result, optimum)
    if result is not None:
        assert len(set(result)) == len(result)
        assert sum(1 << v for v in result) in robust


if __name__ == '__main__':
    count = 0
    for n in range(1, 5):
        options = [mask for mask in range(1 << n) if mask.bit_count() <= 2]
        for rows in product(options, repeat=n):
            check(rows)
            count += 1
    print('Exhaustive looped outdegree <= 2 graphs through n=4:', count, flush=True)
    rng = Random(20261104)
    for n in range(5, 10):
        for _ in range(100):
            rows = []
            for v in range(n):
                neighbors = rng.sample(range(n), rng.randrange(3))
                rows.append(sum(1 << w for w in neighbors))
            check(tuple(rows))
    # Quota overlap zero, several source SCCs, and a rank-three bottleneck.
    fixtures = [(6, 1, 24, 4, 4), (6, 1, 8, 8)]
    fixtures.append(tuple(fixtures[1]) + tuple(r << 4 for r in (6, 1, 1)))
    rank_three = [0] * 9
    for i in range(3):
        c, a, b = 3 * i, 3 * i + 1, 3 * i + 2
        rank_three[c] = (1 << a) | (1 << b)
        rank_three[a] = (1 << c) | (1 << (3 * ((i + 1) % 3)))
        rank_three[b] = 1 << c
    fixtures.append(tuple(rank_three))
    for rows in fixtures:
        check(rows)
    for bad in ((), (-1,), (2,), (7, 0, 0)):
        try:
            solve(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(('invalid input accepted', bad))
    print('Random n=5,...,9 graphs: 500; quota fixtures: 4; input checks: PASS; seed: 20261104')
