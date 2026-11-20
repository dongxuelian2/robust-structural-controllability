"""Exhaustive comparison of SMAWK transitions with direct row scans."""
from itertools import product, combinations_with_replacement

from truncated_partition_fast import _leftmost_row_minima, _vulnerable_transition

checks = 0
for Q in range(1, 6):
    for m in range(Q + 1):
        for slopes in combinations_with_replacement(range(3), m):
            prefix = [0]
            for w in slopes:
                prefix.append(prefix[-1] + w)
            for dp in product(range(4), repeat=Q + 1):
                reference = [
                    min(
                        (dp[u] + prefix[z-u], u)
                        for u in range(max(0, z-m), z+1)
                    )
                    for z in range(Q + 1)
                ]
                penalty = max(dp) + prefix[-1] + 1
                assert _vulnerable_transition(dp, prefix, penalty) == reference
                checks += 1

penalty = 10
transition = _vulnerable_transition([0, penalty, penalty, penalty], [0, 2], penalty)
assert transition[:2] == [(0, 0), (2, 0)]
assert [cost for cost, _ in transition[2:]] == [penalty, penalty]

evaluations = 0


def squared_distance(row, column):
    global evaluations
    evaluations += 1
    return (row - column) ** 2


size = 4096
assert _leftmost_row_minima(size, squared_distance) == list(range(size))
assert evaluations <= 20 * size
print("PASS", checks, "complete small SMAWK min-plus transition checks")
print("PASS unreachable states and", evaluations, "evaluations on a", size, "square implicit matrix")
