"""Independent audit of weighted_fault_core; Python 3.10+, no third parties.

Checks 1,200 seeded binary-matroid instances using exhaustive subset costs
and a worst-case-deletion DP, plus an exact 85-element threshold example.
Run: python3 audit_weighted_independent.py
"""
from functools import lru_cache
from fractions import Fraction
from itertools import combinations
from math import comb
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
    assert set(range(6)) <= core and len(core) == 15
    assert rank_three_one_failure_optimum(labels, weights, range(85)) == 3
    assert rank_three_one_failure_optimum(labels, weights, core) == 3
    # Demonstrate that discarding rejected cheap points at saturation fails.
    uniform, seen = [], set()
    for i in sorted(range(85), key=lambda i: (weights[i], i)):
        if labels[i] not in seen:
            seen.add(labels[i])
            uniform.append(i)
        if len(uniform) == 12:
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


def check_preserver_distinction():
    """Preserving cheapest surviving ordinary bases need not preserve robust OPT.

    Rank-two truncation with parallel classes {0,1}, {2,3}, {4} and
    costs 2,2,2,2,3. The first four points preserve every minimum basis
    under zero or one ground-element failure, but lose a cheaper robust set.
    """
    labels = (0, 0, 1, 1, 2)
    weights = (2, 2, 2, 2, 3)
    ground, preserver = frozenset(range(5)), frozenset(range(4))

    def rank(elements):
        return min(2, len({labels[e] for e in elements}))

    def minimum_basis_cost(available):
        return min(sum(weights[e] for e in pair)
                   for pair in combinations(sorted(available), 2)
                   if rank(pair) == 2)

    def robust_optimum(available):
        feasible_costs = []
        for size in range(len(available) + 1):
            for candidate in combinations(sorted(available), size):
                selected = frozenset(candidate)
                if rank(selected) == 2 and all(
                        rank(selected - {e}) == 2 for e in selected):
                    feasible_costs.append(sum(weights[e] for e in selected))
        return min(feasible_costs)

    failures = [frozenset()] + [frozenset({e}) for e in ground]
    for failure in failures:
        assert minimum_basis_cost(ground - failure) == 4
        assert minimum_basis_cost(preserver - failure) == 4
    assert robust_optimum(ground) == 7
    assert robust_optimum(preserver) == 8


def check_hyperplane_bound():
    """Audit the new count directly, without calling the core or solver.

    Uniform matroids attain equality, including k=0 and rank one.
    On the nonrepresentable Vamos matroid and the Fano matroid, inspect
    every Y and all deletions of at most two elements. This tests the
    more general counting lemma, without assuming Y=B-{x} or resilience.
    """
    def inspect(rank, n, selected, k):
        ground = frozenset(range(n))
        r = rank(ground)
        planes = set()
        for size in range(min(k, len(selected)) + 1):
            for failed in combinations(sorted(selected), size):
                survivor = selected - frozenset(failed)
                if rank(survivor) == r - 1:
                    planes.add(frozenset(e for e in ground
                                        if rank(survivor | {e}) == r - 1))
        pairs = []
        for plane in planes:
            basis = frozenset()
            for e in sorted(selected & plane):
                if rank(basis | {e}) > rank(basis):
                    basis |= {e}
            deleted = selected - plane
            assert len(basis) == r - 1 and len(deleted) <= k
            assert not basis & deleted
            pairs.append((basis, deleted))
        for i, (basis, _) in enumerate(pairs):
            assert all(basis & deleted for j, (_, deleted) in enumerate(pairs)
                       if i != j)
        probability_sum = sum((Fraction(1, comb(r - 1 + len(deleted), r - 1))
                               for _, deleted in pairs), Fraction(0))
        assert probability_sum <= 1
        assert len(planes) <= comb(r + k - 1, k)
        return len(planes)

    for r in range(1, 6):
        for k in range(4):
            n = r + k
            rank = lambda s, r=r: min(r, len(s))
            assert inspect(rank, n, frozenset(range(n - 1)), k) == comb(r+k-1, k)

    forbidden = {frozenset(s) for s in (
        (0,1,2,3), (0,1,4,5), (0,1,6,7), (2,3,4,5), (2,3,6,7))}
    def vamos_rank(s):
        return 3 if s in forbidden else min(4, len(s))
    fano_table = binary_ranks(tuple(range(1, 8)))
    def fano_rank(s):
        return fano_table[sum(1 << e for e in s)]
    for n, rank in ((8, vamos_rank), (7, fano_rank)):
        for mask in range(1 << n):
            selected = frozenset(e for e in range(n) if mask >> e & 1)
            for k in range(3):
                inspect(rank, n, selected, k)



def check_rank_one_shortcut():
    # Compare against the independent deletion-DP on loops, zero weights,
    # safe/vulnerable mixtures, ties and infeasible failure budgets.
    for columns in ((1, 1, 1, 1), (0, 1, 0, 1)):
        for weights in ((0, 0, 0, 0), (1, 2, 3, 4), (4, 3, 2, 1)):
            for k in range(5):
                for mask in range(16):
                    vulnerable = {e for e in range(4) if mask >> e & 1}
                    check_small(columns, weights, k, vulnerable)
    # Enumeration previously visited 2**80 candidates before finding this
    # unique feasible set. Keep this a regression for the direct formula.
    calls = 0
    def rank(s):
        nonlocal calls
        calls += 1
        return int(bool(s - {80}))  # element 80 is a zero-cost loop
    weights = [1] * 80 + [0]
    result = solve(weights, rank, 79)
    assert result[:2] == (80, frozenset(range(80)))
    assert calls <= 100
    # A safe nonloop can win or lose against the vulnerable alternative.
    for safe_cost, expected in ((79, 79), (81, 80)):
        result = solve([1] * 80 + [safe_cost], lambda s: int(bool(s)),
                       79, range(80))
        assert result[0] == expected
    assert solve([1] * 80, lambda s: int(bool(s)), 80)[0] is None


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
    check_rank_one_shortcut()
    check_preserver_distinction()
    check_hyperplane_bound()
    print(f'PASS: {count} seeded instances; {feasible} feasible; '
          f'{reduced} strictly reduced cores.')
    print(f'PASS: 85-element threshold instance; core={core_size}; '
          'original/core optimum=3; uniform-only optimum=13.')
    print('PASS: rank-one shortcut, 480 exhaustive fixtures and large failure budget.')
    print('PASS: iterator inputs, rank zero, infeasibility, and validation.')
    print('PASS: minimum-basis preserver costs remain 4; robust optima differ (7 vs 8).')
    print('PASS: sharp set-pairs hyperplane bound, uniform equality, and '
          'all selected subsets of Vamos/Fano fixtures.')


if __name__ == '__main__':
    main()

