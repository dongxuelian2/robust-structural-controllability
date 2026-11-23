"""Finite exact audits of coefficient compression and one-sided graph sampling."""

from fractions import Fraction
from itertools import combinations, product
from random import Random

from linear_accessibility_quotas import solve_represented
from randomized_actuator_placement import sampled_dual, sampling_bits, solve, verify_placement


def subsets(elements, limit=None):
    elements = tuple(elements)
    maximum = len(elements) if limit is None else min(limit, len(elements))
    for size in range(maximum + 1):
        yield from map(frozenset, combinations(elements, size))


def small_rank(columns):
    """Rank by reducing vectors against one nonzero pivot coordinate at a time."""
    vectors = [tuple(map(Fraction, vector)) for vector in columns if any(vector)]
    if not vectors:
        return 0
    head, *tail = vectors
    j = next(i for i, value in enumerate(head) if value)
    return 1 + small_rank(tuple(a - b * vector[j] / head[j] for a, b in zip(vector, head))
                          for vector in tail)


def optimum(n, weights, k, vulnerable, survives):
    best = None
    for chosen in subsets(range(n)):
        cost = sum(weights[e] for e in chosen)
        if best is None or cost < best:
            if all(survives(chosen - failure) for failure in subsets(chosen & vulnerable, k)):
                best = cost
    return best


def represented_audit():
    rng = Random(1121)
    checked = 0
    for n in range(1, 7):
        for _ in range(12):
            dimension = rng.randrange(3)
            columns = tuple(tuple(rng.randrange(-2, 3) for _ in range(dimension))
                            for e in range(n))
            weights = tuple(rng.randrange(5) for e in range(n))
            vulnerable = frozenset(e for e in range(n) if rng.randrange(2))
            colors = [rng.randrange(3) for e in range(n)]
            groups = [frozenset(e for e in range(n) if colors[e] == g)
                      for g in sorted(set(colors) - {0})]
            k = rng.randrange(3)
            r = small_rank(columns)

            def survives(survivor):
                return all(survivor & group for group in groups) and small_rank(
                    columns[e] for e in survivor
                ) == r

            expected = optimum(n, weights, k, vulnerable, survives)
            actual = solve_represented(columns, weights, k, vulnerable, groups)
            assert (None if actual is None else sum(weights[e] for e in actual.placement)) == expected
            if actual is not None:
                assert all(survives(actual.placement - failure)
                           for failure in subsets(actual.placement & vulnerable, k))
                assert actual.witness <= actual.placement and len(actual.witness) <= (k + 1) * r
                assert all(small_rank(columns[e] for e in actual.witness - failure) == r
                           for failure in subsets(actual.witness & vulnerable, k))
            checked += 1

    # The three singleton groups force quota-only loops beyond the witness bound.
    columns = [(Fraction(2, 3), 0), (0, 0), (0, 0), (0, 0), (3, 0)]
    groups = [{1}, {2}, {3}]
    solution = solve_represented(columns, (2, 3, 4, 5, 1), 1, {4}, groups)
    assert solution.placement == {0, 1, 2, 3} and solution.witness == {0}
    assert solve_represented([(0,)], (0,), 0, groups=[set()]) is None

    columns = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (1, -1, 2), (0, 0, 0)]
    weights = (1, 2, 4, 2, 1, 1)
    exposed = frozenset({0, 1, 3, 4})
    groups = [{0, 1, 5}, {2, 3, 4}]

    def survives(survivor):
        return all(survivor & group for group in groups) and small_rank(
            columns[e] for e in survivor
        ) == 3

    expected = optimum(len(columns), weights, 1, exposed, survives)
    solution = solve_represented(columns, weights, 1, exposed, groups)
    assert sum(weights[e] for e in solution.placement) == expected
    assert all(survives(solution.placement - failure)
               for failure in subsets(solution.placement & exposed, 1))
    checked += 1
    return checked


def can_match(rows, targets):
    if not targets:
        return True
    target, *rest = targets
    return any(row >> target & 1 and can_match(rows[:i] + rows[i + 1:], rest)
               for i, row in enumerate(rows))


def graph_survives(rows, selected):
    n = len(rows)
    reached = set(selected)
    front = list(selected)
    for v in front:
        for u in range(n):
            if rows[v] >> u & 1 and u not in reached:
                reached.add(u)
                front.append(u)
    return len(reached) == n and can_match(rows, sorted(set(range(n)) - selected))


def basis_inclusion_audit():
    accepted = drops = 0
    for masks in product(range(4), repeat=2):
        edges = [(i, j) for i in range(2) for j in range(2) if masks[i] >> j & 1]
        primal = [choice for choice in subsets(range(2)) if can_match(masks, sorted(choice))]
        rank = max(map(len, primal))
        for entries in product((1, 2), repeat=len(edges)):
            matrix = [[0, 0], [0, 0]]
            for (i, j), value in zip(edges, entries):
                matrix[i][j] = value
            columns = sampled_dual(matrix, rank)
            if columns is None:
                assert small_rank(zip(*matrix)) < rank
                drops += 1
                continue
            accepted += 1
            assert small_rank(columns) == 2 - rank
            for choice in subsets(range(2)):
                if small_rank(columns[e] for e in choice) == 2 - rank:
                    assert frozenset(range(2)) - choice in primal
    return accepted, drops


class _FixedDraws:
    def __init__(self, values):
        self.values = iter(values)

    def getrandbits(self, bits):
        value = next(self.values) - 1
        assert 0 <= value < 1 << bits
        return value


def graph_audit():
    rng = Random(1122)
    checked = 0
    for n in range(1, 5):
        for _ in range(18):
            rows = tuple(rng.randrange(1 << n) for i in range(n))
            vulnerable = frozenset(e for e in range(n) if rng.randrange(2))
            weights = tuple(rng.randrange(6) for e in range(n))
            k = rng.randrange(3)
            exact = optimum(n, weights, k, vulnerable, lambda chosen: graph_survives(rows, chosen))
            trial = solve(rows, weights, k, vulnerable, error_bits=30, rng=rng)
            cost = None if trial.placement is None else sum(weights[e] for e in trial.placement)
            assert cost == exact
            if trial.placement is not None:
                assert verify_placement(rows, k, vulnerable, trial.placement, trial.witness)
                assert all(graph_survives(rows, trial.placement - failure)
                           for failure in subsets(trial.placement & vulnerable, k))
            checked += 1

    rows = (24, 24, 24, 31, 31)
    matrix = [[0, 0, 0, 1, 1], [0, 0, 0, 1, 2], [0, 0, 0, 2, 1],
              [1, 1, 1, 1, 1], [1, 1, 2, 1, 1]]
    weights = (4, 6, 1, 9, 9)
    exposed = frozenset({0, 1, 3, 4})
    assert graph_survives(rows, {2})
    for k, expected in ((1, {0, 1}), (2, None)):
        draws = _FixedDraws(value for row in matrix for value in row if value)
        trial = solve(rows, weights, k, exposed, rng=draws)
        assert trial.placement == expected
        assert trial.status == ("candidate" if expected else "no_candidate")
        if expected:
            assert verify_placement(rows, k, exposed, trial.placement, trial.witness)
            assert not verify_placement(rows, k, exposed, trial.placement, {0})
    draws = _FixedDraws(value for row in matrix for value in row if value)
    equal = solve(rows, (1, 4, 9, 9, 9), 1, set(), rng=draws)
    assert equal.placement == {0}
    draws = _FixedDraws([1] * sum(row.bit_count() for row in rows))
    assert solve(rows, weights, 1, exposed, rng=draws).status == "rank_drop"

    assert solve((0,), (1,), 1, {0}, rng=_FixedDraws([])).placement is None
    assert solve((0,), (1,), 1, set(), rng=_FixedDraws([])).placement == {0}
    assert solve((1,), (1,), 1, set(), rng=_FixedDraws([])).placement == {0}
    assert sampling_bits(1000, 998, 2, 1000, 20) == 33
    return checked


if __name__ == "__main__":
    print(f"PASS: {represented_audit()} represented-instance optimum comparisons.")
    accepted, drops = basis_inclusion_audit()
    print(f"PASS: {accepted} rank-preserving samples and {drops} rank drops checked.")
    print(f"PASS: {graph_audit()} graph optima, rank-preserving misses, and rank-zero cases.")
