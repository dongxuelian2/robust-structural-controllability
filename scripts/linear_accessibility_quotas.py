"""Exact rational coefficient DP for weighted resilient spanning and group survival."""

from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, permutations
from math import lcm
from typing import NamedTuple

from weighted_fault_core import weighted_input


def reduced_rows(rows):
    """Return the reduced nonzero rows and their pivot columns over Q."""
    matrix = [list(map(Fraction, row)) for row in rows]
    if not matrix:
        return [], ()
    pivots = []
    for column in range(len(matrix[0])):
        index = len(pivots)
        pivot = next((i for i in range(index, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[index], matrix[pivot] = matrix[pivot], matrix[index]
        divisor = matrix[index][column]
        matrix[index] = [value / divisor for value in matrix[index]]
        for i, row in enumerate(matrix):
            if i != index and row[column]:
                factor = row[column]
                matrix[i] = [a - factor * b for a, b in zip(row, matrix[index])]
        pivots.append(column)
        if len(pivots) == len(matrix):
            break
    return matrix[:len(pivots)], tuple(pivots)


def _failures(elements, k):
    for size in range(min(len(elements), k) + 1):
        yield from combinations(elements, size)


def _add_scaled(result, polynomial, scale=1):
    for term, coefficient in polynomial.items():
        value = result.get(term, 0) + scale * coefficient
        if value:
            result[term] = value
        else:
            result.pop(term, None)


def _multiply(left, right):
    # A monomial is a sorted tuple of variable indices with repetition.
    result = defaultdict(int)
    for x, a in left.items():
        for y, b in right.items():
            result[tuple(sorted(x + y))] += a * b
    return {term: value for term, value in result.items() if value}


@cache
def _signed_permutations(size):
    return tuple(
        (order, (-1) ** sum(a > b for i, a in enumerate(order) for b in order[i + 1:]))
        for order in permutations(range(size))
    )


def _determinant(columns):
    result = {}
    for order, sign in _signed_permutations(len(columns)):
        term = {(): sign}
        for column, row in zip(columns, order):
            term = _multiply(term, column[row])
            if not term:
                break
        _add_scaled(result, term)
    return result


def _certificate(witness, completed_counts, columns, basis, vulnerable, rank, k):
    """Common-completion polynomial with distinct integer tags for future slots."""
    future_safe, future_vulnerable = completed_counts
    slots = future_safe + future_vulnerable
    vulnerable_slots = range(future_safe, slots)
    exposed = tuple(e for e in witness if e in vulnerable)
    product = {(): 1}
    for past_failure in _failures(exposed, k):
        past = basis(tuple(e for e in witness if e not in past_failure))
        missing = rank - len(past)
        if missing == 0:
            continue
        fixed = [[{(): value} if value else {} for value in columns[e]] for e in past]
        for future_failure in _failures(vulnerable_slots, k - len(past_failure)):
            if missing > slots - len(future_failure):
                return {}
            coefficients = [1]
            for slot in future_failure:
                tag = slot + 1
                coefficients = (
                    [-tag * coefficients[0]]
                    + [a - tag * b for a, b in zip(coefficients, coefficients[1:])]
                    + [coefficients[-1]]
                )
            moments = [
                [{((h + j) * rank + row,): c for j, c in enumerate(coefficients) if c}
                 for row in range(rank)]
                for h in range(missing)
            ]
            product = _multiply(product, _determinant(fixed + moments))
            if not product:
                return {}
    return product


class _Purchase(NamedTuple):
    cost: int
    selected: tuple
    witness: tuple


class PurchaseSolution(NamedTuple):
    placement: frozenset
    witness: frozenset


def _compress(candidates, vector):
    pivots = {}
    kept = []
    for candidate in sorted(candidates):
        row = vector(candidate.witness).copy()
        while row:
            first = min(row)
            if first in pivots:
                _add_scaled(row, pivots[first], -row[first])
            else:
                pivot = Fraction(row[first])
                pivots[first] = {term: value / pivot for term, value in row.items()}
                kept.append(candidate)
                break
    return kept


def _cheapest_group_cover(groups, weights, vulnerable, threshold):
    selected = set()
    for group in groups:
        safe = sorted(group - vulnerable, key=lambda e: (weights[e], e))
        exposed = sorted(group & vulnerable, key=lambda e: (weights[e], e))
        options = [safe[:1]] if safe else []
        if len(exposed) >= threshold:
            options.append(exposed[:threshold])
        if not options:
            return None
        selected.update(min(options, key=lambda seq: sum(weights[e] for e in seq)))
    return frozenset(selected)


def solve_represented(columns, weights, k, vulnerable=None, groups=()):
    """Return a cheapest feasible purchase and its small witness, or None.

    Elements are indexed by rational column. Each group must retain a purchase
    after any failure of at most k selected vulnerable elements.
    """
    weights, ground, vulnerable = weighted_input(weights, k, vulnerable)
    columns = tuple(tuple(map(Fraction, column)) for column in columns)
    if len(columns) != len(weights) or len({len(col) for col in columns}) > 1:
        raise ValueError("Expected one column of the same dimension per element.")
    groups = tuple(map(frozenset, groups))
    grouped = set()
    for group in groups:
        if grouped & group or not group <= ground:
            raise ValueError("Groups must be disjoint subsets of the ground set.")
        grouped.update(group)
    k = min(k, len(vulnerable))
    threshold = k + 1
    group_cover = _cheapest_group_cover(groups, weights, vulnerable, threshold)
    if group_cover is None:
        return None

    row_basis, _ = reduced_rows(zip(*columns))
    rank = len(row_basis)
    if rank == 0:
        return PurchaseSolution(group_cover, frozenset())
    # Scaling each column by a nonzero denominator preserves its matroid.
    integral = []
    for column in zip(*row_basis):
        scale = lcm(*(value.denominator for value in column))
        integral.append(tuple(int(value * scale) for value in column))
    columns = tuple(integral)

    @cache
    def basis(elements):
        _, pivots = reduced_rows(zip(*(columns[e] for e in elements)))
        return tuple(elements[i] for i in pivots)

    def solve_target(total_safe, total_vulnerable):
        @cache
        def certificate(witness, sa, va):
            return _certificate(witness, (total_safe - sa, total_vulnerable - va),
                                columns, basis, vulnerable, rank, k)

        def compress_buckets(buckets):
            return {
                state: kept
                for state, pool in buckets.items()
                if (kept := _compress(pool, lambda w: certificate(w, state[0], state[1])))
            }

        states = {(0, 0, 0, False, 0): [_Purchase(0, (), ())]}
        for position, block in enumerate((*groups, ground - grouped)):
            requires_cover = position < len(groups)
            for e in sorted(block):
                exposed = e in vulnerable
                buckets = defaultdict(list)
                for state, pool in states.items():
                    sa, va, p, has_safe, bought_vulnerable = state
                    buckets[state].extend(pool)
                    new_safe = has_safe or (requires_cover and not exposed)
                    new_vulnerable = min(threshold, bought_vulnerable + (requires_cover and exposed))
                    for candidate in pool:
                        purchases = candidate.selected + (e,)
                        cost = candidate.cost + weights[e]
                        buckets[sa, va, p, new_safe, new_vulnerable].append(
                            _Purchase(cost, purchases, candidate.witness)
                        )
                        if (exposed and va < total_vulnerable) or (not exposed and sa < total_safe):
                            witness = candidate.witness + (e,)
                            next_state = (sa + (not exposed), va + exposed, len(basis(witness)),
                                          new_safe, new_vulnerable)
                            buckets[next_state].append(_Purchase(cost, purchases, witness))
                states = compress_buckets(buckets)
            if requires_cover:
                passed = defaultdict(list)
                for (sa, va, p, has_safe, bought_vulnerable), pool in states.items():
                    if has_safe or bought_vulnerable == threshold:
                        passed[sa, va, p, False, 0].extend(pool)
                states = compress_buckets(passed)
        return min(states.get((total_safe, total_vulnerable, rank, False, 0), ()), default=None)

    limit = min(len(ground), threshold * rank)
    best = None
    for sa in range(min(len(ground - vulnerable), limit) + 1):
        for va in range(max(0, rank - sa), min(len(vulnerable), limit - sa) + 1):
            candidate = solve_target(sa, va)
            if candidate is not None and (best is None or candidate < best):
                best = candidate
    if best is None:
        return None
    return PurchaseSolution(frozenset(best.selected), frozenset(best.witness))
