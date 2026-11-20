"""Exhaustive selected-set/failure test for oracle direct-sum recovery."""
from itertools import combinations
from random import Random
from direct_sum_partition_fault import solve

rng = Random(2026111714)
cases = 0
for t in range(700):
    blocks = rng.randrange(1, 5)
    cfg = []
    nums = []
    comps = []
    for blk in range(blocks):
        p = rng.randrange(1, 4)
        sizes = [rng.randrange(1, 3) for _ in range(p)]
        rr = rng.randrange(p+1)
        for cls, size in enumerate(sizes):
            for _ in range(size):
                comps.append(blk)
                nums.append(cls)
        cfg.append(rr)
    if len(nums) > 11:
        continue
    for _ in range(rng.randrange(2)):
        comps.append(blocks)
        nums.append(-1)
    n = len(nums)
    weights = [rng.randrange(9) for _ in range(n)]
    vulnerable = frozenset(i for i in range(n) if rng.randrange(2))
    k = rng.choice([0, 1, 2, 3, n+2])

    def rank(s):
        ans = 0
        for blk, rr in enumerate(cfg):
            ans += min(rr, len(set(nums[x] for x in s
                                   if comps[x] == blk and nums[x] >= 0)))
        return ans

    R = sum(cfg)

    def independent(s):
        return len(s) == rank(s)

    val, witness, info = solve(weights, independent, k, vulnerable)
    brute_optimum = None
    for mask in range(1 << n):
        B = frozenset(i for i in range(n) if mask >> i & 1)
        cost = sum(weights[x] for x in B)
        if brute_optimum is not None and cost >= brute_optimum:
            continue
        can_fail = sorted(B & vulnerable)
        if all(rank(B-set(F)) == R
               for f in range(min(k, len(can_fail)) + 1)
               for F in combinations(can_fail, f)):
            brute_optimum = cost
    assert val == brute_optimum, (
        t, n, cfg, nums, comps, weights, k, vulnerable,
        val, brute_optimum, witness, info
    )
    if witness is not None:
        assert sum(weights[x] for x in witness) == val
        assert all(rank(witness-set(F)) == R
                   for f in range(min(k, len(witness & vulnerable)) + 1)
                   for F in combinations(sorted(witness & vulnerable), f))
    cases += 1

print("PASS", cases, "direct-sum oracle instances against exhaustive optima")
