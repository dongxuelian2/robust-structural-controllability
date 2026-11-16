"""Exhaustive finite Monge-row-minima check for arbitrary predecessor arrays.

This checks the staircase interval bounds, monotonicity of the leftmost
argmin, and divide-and-conquer minima on every tested array/convex prefix.
It is independent of the matroid solver and not a mathematical proof.
"""
from itertools import product, combinations_with_replacement

checks = 0
for Q in range(1, 6):
    for m in range(Q + 1):
        for slopes in combinations_with_replacement(range(3), m):
            prefix = [0]
            for w in slopes:
                prefix.append(prefix[-1] + w)
            for dp in product(range(4), repeat=Q + 1):
                reference = []
                argmins = []
                for z in range(Q + 1):
                    val, arg = min(
                        (dp[u] + prefix[z-u], u)
                        for u in range(max(0, z-m), z+1)
                    )
                    reference.append(val)
                    argmins.append(arg)
                assert argmins == sorted(argmins), (Q, m, prefix, dp, argmins)

                accelerated = [None] * (Q + 1)

                def divide_and_conquer(lo, hi, left, right):
                    if lo > hi:
                        return
                    z = (lo + hi) // 2
                    val, arg = min(
                        (dp[u] + prefix[z-u], u)
                        for u in range(max(left, z-m), min(right, z)+1)
                    )
                    accelerated[z] = val
                    divide_and_conquer(lo, z-1, left, arg)
                    divide_and_conquer(z+1, hi, arg, right)

                divide_and_conquer(0, Q, 0, Q)
                assert accelerated == reference
                checks += 1

print("PASS", checks, "complete small Monge min-plus transition checks")
