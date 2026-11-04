"""Independent audit of weighted_fault_core; Python 3.10+, no third parties.

Checks 1,200 seeded binary-matroid instances using exhaustive subset costs
and a worst-case-deletion DP, plus an exact 85-element threshold example.
Run: python3 audit_weighted_independent.py
"""
from functools import lru_cache
from random import Random

from weighted_fault_core import solve, weighted_core


def binary_ranks(columns):
    """Independent GF(2) elimination for every subset of the column list."""
    ranks = []
    for mask in range(1 << len(columns)):
        basis = {}
        for i, column in enumerate(columns):
            if not mask >> i & 1:
                continue
            value = column
            while value:
                pivot = value.bit_length() - 1
                if pivot in basis:
                    value ^= basis[pivot]
                else:
                    basis[pivot] = value
                    break
        ranks.append(len(basis))
    return ranks


def check_small(columns, weights, k, vulnerable):
    ranks = binary_ranks(columns)
    n = len(weights)
    full = (1 << n) - 1
    target = ranks[full]
    vulnerable_mask = sum(1 << i for i in vulnerable)
    # After j iterations, survival[s] is the minimum rank obtainable by
    # deleting at most j vulnerable elements from s. This differs from the
    # candidate implementation's explicit fixed-cardinality failure checks.
    survival = ranks[:]
    for _ in range(k):
        previous = survival
        survival = [
            min([previous[s]] + [previous[s ^ (1 << i)]
                 for i in range(n) if s >> i & 1 and vulnerable_mask >> i & 1])
            for s in range(full + 1)
        ]
    costs = [0] * (full + 1)
    for s in range(1, full + 1):
        bit = s & -s
        costs[s] = costs[s ^ bit] + weights[bit.bit_length() - 1]
    optimum = min((costs[s] for s in range(full + 1)
                   if survival[s] == target), default=None)

    def oracle(elements):
        return ranks[sum(1 << i for i in elements)]

    actual, witness, core, _ = solve(iter(weights), oracle, k, iter(vulnerable))
    assert actual == optimum, (columns, weights, k, vulnerable, actual, optimum)
    if witness is not None:
        mask = sum(1 << i for i in witness)
        assert witness <= core
        assert costs[mask] == actual and survival[mask] == target
    else:
        assert actual is None
    return len(core) < n, actual is not None


class TruncatedParallelRanks:
    """Rank-r truncation of a partition matroid with unit class capacities."""
    def __init__(self, labels, rank):
        self.labels = tuple(labels)
        self.rank = rank

    @lru_cache(None)
    def __call__(self, elements):
        return min(self.rank, len({self.labels[i] for i in elements}))


def rank_three_one_failure_optimum(labels, weights, allowed):
    """Exact nonnegative-cost optimum for a rank-three parallel truncation.

    Any feasible set with >=4 occupied classes contains a feasible choice
    of four singleton representatives. With exactly three occupied classes,
    each must contain >=2 elements. Fewer classes cannot span robustly.
    """
    classes = {}
    for i in allowed:
        classes.setdefault(labels[i], []).append(weights[i])
    for costs in classes.values():
        costs.sort()
    singles = sorted(costs[0] for costs in classes.values())
    doubles = sorted(sum(costs[:2]) for costs in classes.values() if len(costs) >= 2)
    candidates = []
    if len(singles) >= 4:
        candidates.append(sum(singles[:4]))
    if len(doubles) >= 3:
        candidates.append(sum(doubles[:3]))
    return min(candidates, default=None)


def check_threshold():
    labels = [0, 0, 1, 1, 2, 2] + list(range(3, 82))
    weights = [0, 1, 0, 1, 0, 1] + [10 + c for c in range(3, 82)]
    oracle = TruncatedParallelRanks(labels, 3)
    core, stats = weighted_core(weights, oracle, 1)
    assert stats['capped'] >= 1 and stats['rejected'] >= 1
    assert set(range(6)) <= core and len(core) == 29
    assert rank_three_one_failure_optimum(labels, weights, range(85)) == 3
    assert rank_three_one_failure_optimum(labels, weights, core) == 3
    # Demonstrate that discarding rejected cheap points at saturation fails.
    uniform, seen = [], set()
    for i in sorted(range(85), key=lambda i: (weights[i], i)):
        if labels[i] not in seen:
            seen.add(labels[i])
            uniform.append(i)
        if len(uniform) == 26:
            break
    assert rank_three_one_failure_optimum(labels, weights, uniform) == 13
    return len(core)


def check_edges():
    assert solve([], lambda s: 0, 100)[:2] == (0, frozenset())
    assert solve([0], lambda s: len(s), 0)[0] == 0
    assert solve([1], lambda s: len(s), 1)[0] is None
    assert solve([1], lambda s: len(s), 10, vulnerable=[])[0] == 1
    invalid = (([-1], 0, None), ([True], 0, None), ([1], -1, None),
               ([1], False, None), ([1], 0, [2]))
    for weights, k, vulnerable in invalid:
        try:
            weighted_core(weights, lambda s: len(s), k, vulnerable)
        except ValueError:
            pass
        else:
            raise AssertionError(('invalid input accepted', weights, k, vulnerable))


def main():
    rng = Random(202611041236)
    count = reduced = feasible = 0
    for n in range(1, 13):
        for _ in range(100):
            limit = rng.randrange(1, 5)
            columns = [rng.randrange(1 << limit) for _ in range(n)]
            weights = [rng.randrange(31) for _ in range(n)]
            k = rng.randrange(3)
            vulnerable = (set(range(n)) if rng.randrange(2)
                          else {i for i in range(n) if rng.randrange(2)})
            reduction, exists = check_small(columns, weights, k, vulnerable)
            count += 1
            reduced += reduction
            feasible += exists
    core_size = check_threshold()
    check_edges()
    print(f'PASS: {count} seeded instances; {feasible} feasible; '
          f'{reduced} strictly reduced cores.')
    print(f'PASS: 85-element threshold instance; core={core_size}; '
          'original/core optimum=3; uniform-only optimum=13.')
    print('PASS: iterator inputs, rank zero, infeasibility, and validation.')


if __name__ == '__main__':
    main()
