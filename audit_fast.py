"""Cross-check Monge-accelerated solver versus original DP and failures."""
from itertools import combinations
from random import Random
from truncated_partition_fault import solve as slow
from truncated_partition_fast import solve as fast

rng = Random(20480916)
checks = 0
for i in range(1500):
    p = rng.randrange(1, 9)
    sizes = [rng.randrange(1, 4) for _ in range(p)]
    if sum(sizes) > 13:
        continue
    ids = [v for v, size in enumerate(sizes) for _ in range(size)]
    ids += [-1] * rng.randrange(0, 3)
    r = rng.randrange(p+1)
    n = len(ids)
    weights = [rng.randrange(15) for _ in range(n)]
    v = frozenset(j for j in range(n) if rng.randrange(2))
    k = rng.choice([0, 1, 2, 3, 4, n, n+4])

    def independent(s):
        a = [ids[x] for x in s]
        return len(a) <= r and -1 not in a and len(a) == len(set(a))

    def rank(s):
        return min(r, len(set(ids[j] for j in s if ids[j] >= 0)))

    a = slow(weights, independent, k, v)
    b = fast(weights, independent, k, v)
    assert a[0] == b[0], (i, ids, r, k, weights, v, a, b)
    if b[1] is not None:
        assert sum(weights[j] for j in b[1]) == b[0]
        assert all(rank(b[1] - set(F)) == r
                   for d in range(min(k, len(b[1] & v)) + 1)
                   for F in combinations(sorted(b[1] & v), d))
    checks += 1

print("PASS", checks, "seeded solver comparisons and failure certificates")
