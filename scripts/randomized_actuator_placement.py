"""One-sided rational sampling for weighted dedicated actuator placement."""

from itertools import combinations
from math import comb
from random import SystemRandom
from typing import NamedTuple

from bounded_core_solver import matching_number, source_components
from linear_accessibility_quotas import reduced_rows, solve_represented
from weighted_fault_core import weighted_input


class Trial(NamedTuple):
    status: str
    placement: frozenset | None
    witness: frozenset | None
    sample_bits: int


def sampling_bits(n, primal_rank, k, vulnerable_count, error_bits):
    """Size the sampling interval for optimality with error at most 2**-error_bits."""
    if primal_rank == 0 or primal_rank == n:
        return 0
    witness_size = min(n, (k + 1) * (n - primal_rank))
    exposed = min(witness_size, vulnerable_count)
    failure_cap = min(k, exposed)
    failures = sum(comb(exposed, j) for j in range(failure_cap + 1))
    basis_cover = comb(n - primal_rank + failure_cap, failure_cap)
    minors = min(comb(n, primal_rank), failures, basis_cover)
    return error_bits + (primal_rank * minors - 1).bit_length()


def sampled_dual(matrix, primal_rank):
    """Return rational nullspace columns if the sample has the prescribed rank."""
    reduced, pivots = reduced_rows(matrix)
    if len(pivots) != primal_rank:
        return None
    n = len(matrix[0]) if matrix else 0
    rows = []
    for free in sorted(set(range(n)) - set(pivots)):
        row = [0] * n
        row[free] = 1
        for i, pivot in enumerate(pivots):
            row[pivot] = -reduced[i][free]
        rows.append(row)
    return tuple(zip(*rows)) if rows else ((),) * n


def verify_placement(rows, k, vulnerable, placement, witness):
    """Check accessibility and the witness by state matching on the original graph."""
    rows = tuple(rows)
    ground = frozenset(range(len(rows)))
    vulnerable = ground if vulnerable is None else frozenset(vulnerable)
    if type(k) is not int or k < 0 or placement is None or witness is None:
        return False
    if not witness <= placement <= ground or not vulnerable <= ground:
        return False
    selected_mask = sum(1 << e for e in placement)
    vulnerable_mask = sum(1 << e for e in vulnerable)
    for source in source_components(rows):
        bought = selected_mask & source
        if not bought & ~vulnerable_mask and bought.bit_count() <= k:
            return False
    full = (1 << len(rows)) - 1
    outside = full ^ sum(1 << e for e in witness)
    exposed = sorted(witness & vulnerable)
    for size in range(min(k, len(exposed)) + 1):
        for failed in combinations(exposed, size):
            targets = outside | sum(1 << e for e in failed)
            if matching_number(rows, targets) != targets.bit_count():
                return False
    return True


def solve(rows, weights, k, vulnerable=None, *, error_bits=20, rng=None):
    """Run one trial on a directed state graph in outgoing-bitmask format.

    Status is candidate, rank_drop, or no_candidate. Every returned placement
    is feasible on the input graph. A feasible input yields an optimum with
    probability at least 1 - 2**-error_bits under independent uniform draws.
    """
    weights, ground, vulnerable = weighted_input(weights, k, vulnerable)
    rows = tuple(rows)
    n = len(rows)
    if n != len(weights) or n == 0 or any(type(row) is not int or not 0 <= row < 1 << n
                                          for row in rows):
        raise ValueError("Expected a nonempty graph and one weight per state.")
    if type(error_bits) is not int or error_bits < 1:
        raise ValueError("error_bits must be a positive integer.")
    k = min(k, len(vulnerable))
    primal_rank = matching_number(rows, (1 << n) - 1)
    groups = tuple(frozenset(e for e in ground if mask >> e & 1)
                   for mask in source_components(rows))
    bits = sampling_bits(n, primal_rank, k, len(vulnerable), error_bits)
    if primal_rank == 0:
        if k:
            return Trial("no_candidate", None, None, 0)
        return Trial("candidate", ground, ground, 0)
    if primal_rank == n:
        columns = ((),) * n
    else:
        generator = SystemRandom() if rng is None else rng
        matrix = [[1 + generator.getrandbits(bits) if row >> e & 1 else 0
                   for e in range(n)] for row in rows]
        columns = sampled_dual(matrix, primal_rank)
        if columns is None:
            return Trial("rank_drop", None, None, bits)
    solution = solve_represented(columns, weights, k, vulnerable, groups)
    if solution is None:
        return Trial("no_candidate", None, None, bits)
    return Trial("candidate", solution.placement, solution.witness, bits)
