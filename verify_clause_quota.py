"""Exact clause ranks, global consistency, and sparse control certificates."""
from itertools import combinations, combinations_with_replacement, product
from random import Random
from bounded_core_solver import matching_number, source_components
from clause_quota import CLAUSE_ROWS, clause_quota_instance, lift_truth, round_truth, violated
from quota_cotransversal import quota_cotransversal_graph, quota_from_placement
from verify_gammoid_obstructions import dual_ranks, robust, transversal_ranks
from verify_quota_cotransversal import input_first_matchable
from verify_scc_cover import girth


def clause_masks():
    ranks = dual_ranks(CLAUSE_ROWS, 7)
    assert ranks[-1] == 3
    assert all(ranks[mask] == mask.bit_count() for mask in range(128) if mask.bit_count() <= 2)
    assert {mask for mask in range(128) if mask.bit_count() == 3 and ranks[mask] < 3} == set(CLAUSE_ROWS)
    assert all(ranks[mask] == 3 for mask in range(128) if mask.bit_count() >= 4)
    actual = tuple(mask for mask in range(128) if mask & 1 and robust(ranks, mask))
    patterns = {1 | sum(1 << (2 * i + 1 + bit) for i, bit in enumerate(bits))
                for bits in product((0, 1), repeat=3) if bits != (1, 1, 1)}
    assert set(actual) == patterns | {mask for mask in range(128) if mask & 1 and mask.bit_count() >= 5}
    bases = 0
    for basis in combinations(range(7), 3):
        mask = sum(1 << e for e in basis)
        if ranks[mask] != 3:
            continue
        rows = [(1 << b) | sum(1 << e for e in range(7) if not mask >> e & 1 and
                               ranks[(mask ^ (1 << b)) | (1 << e)] == 3) for b in basis]
        relaxed = transversal_ranks(rows, 7)
        assert any(robust(relaxed, s) and not robust(ranks, s)
                   for s in range(128) if s & 1 and s.bit_count() == 4)
        bases += 1
    assert bases == 31
    print('Seven-element clause: 128 ranks, seven truth patterns, all 31 bases PASS', flush=True)
    return actual


def bitmask(elements):
    return sum(1 << e for e in elements)


def source_valid(instance, selected, local):
    if not instance.forced <= selected or any(len(selected & set(group)) < 2 for group in instance.groups):
        return False
    selected = bitmask(selected)
    return all(((selected >> (7 * c)) & 127) in local for c in range(len(instance.clauses)))


def exact_two_clause(local):
    clauses = tuple(combinations_with_replacement((-2, -1, 1, 2), 3))
    formulas = candidates = feasible = lifts = 0
    for formula in combinations_with_replacement(clauses, 2):
        instance = clause_quota_instance(formula)
        group_masks = tuple(map(bitmask, instance.groups))
        assert sorted(e for group in instance.groups for e in group) == list(range(20))
        assignments = [dict(zip(sorted(instance.ports), bits))
                       for bits in product((False, True), repeat=len(instance.ports))]
        minimum_violations = min(violated(formula, a) for a in assignments)
        for a in assignments:
            lift = lift_truth(instance, a)
            assert source_valid(instance, lift, local)
            lifts += 1
        dummy_mask = bitmask(range(14, 20))
        optimum = 21
        bound = max(map(len, instance.ports.values()))
        loss = 1 + 2 * ((bound - 1) // 2)
        for first, second in product(local, repeat=2):
            chosen = first | (second << 7) | dummy_mask
            candidates += 1
            if any((chosen & group).bit_count() < 2 for group in group_masks):
                continue
            feasible += 1
            support = {e for e in range(20) if chosen >> e & 1}
            assignment = round_truth(instance, support)
            assert violated(formula, assignment) <= loss * (len(support) - 14)
            optimum = min(optimum, len(support))
        assert 14 <= optimum <= 14 + minimum_violations
        assert (optimum == 14) == (minimum_violations == 0)
        assert minimum_violations <= loss * (optimum - 14)
        formulas += 1
    print('Two-clause exact optima:', formulas, 'formulas,', candidates, 'candidate supports,',
          feasible, 'rounded supports,', lifts, 'truth lifts PASS', flush=True)


def cycle_break(k):
    """Three distinct variables per clause, u*=k but excess exactly one."""
    clauses, phase = [], []
    next_var = 2
    for literal in (-1, 1):
        for _ in range(k):
            for a, b in product((1, -1), repeat=2):
                clauses.append((literal, a * next_var, b * (next_var + 1)))
                phase.append(literal > 0)
            next_var += 2
        if literal < 0:
            switch = len(clauses)
            clauses.append((1, next_var, next_var + 1))
            phase.append(False)
            next_var += 2
    clauses.append((next_var, next_var + 1, next_var + 2))
    phase.append(False)
    instance = clause_quota_instance(clauses)
    selected = set(instance.forced)
    for c, clause in enumerate(clauses):
        for j, literal in enumerate(clause):
            value = phase[c] if abs(literal) == 1 else c in (switch, len(clauses) - 1)
            selected.add(7 * c + 2 * j + (1 if value == (literal > 0) else 2))
    selected.add(7 * switch + 1)
    return instance, frozenset(selected)


def check_graph(instance, lifts, local, length):
    graph = quota_cotransversal_graph(instance.presentation, instance.elements, instance.groups, length)
    rows = graph.rows
    n, m = len(rows), len(instance.clauses)
    full, rank = (1 << n) - 1, 9 * m // 2
    assert matching_number(rows, full) == n - rank
    assert graph.core_vertices == 113 * m // 2 - 2
    assert graph.core_edges == 143 * m // 2 - 2
    assert graph.leaves == 0
    sources = source_components(rows)
    assert sorted(s.bit_count() for s in sources) == [4 * length] * (m // 2) + [6 * length] * (3 * m)
    assert 2 * len(sources) == 7 * m
    assert girth(rows) >= 4 * length
    # Each new column comes from one old block. The disjoint target
    # restrictions retain ranks 3, 1 or 0, which add to the whole rank.
    block = [v // 7 for v in range(7 * m)] + [m + v // 2 for v in range(3 * m)]
    for _row, old in graph.steps:
        owners = {block[e] for e in range(len(block)) if old >> e & 1}
        assert len(owners) == 1
        block.append(owners.pop())
    partition = [set() for _ in range(5 * m // 2 + 1)]
    for v, origin in enumerate(graph.projection):
        partition[block[origin] if origin >= 0 else -1].add(v)
    assert [len(part) + matching_number(rows, full ^ bitmask(part)) - (n - rank)
            for part in partition] == [3] * m + [1] * (3 * m // 2) + [0]

    def actual(selected):
        mask = bitmask(selected)
        return all((mask & group).bit_count() >= 2 for group in sources) and all(
            input_first_matchable(rows, selected - {v}) for v in selected)

    for source in lifts:
        assert source_valid(instance, source, local)
        target = {graph.representatives[e] for e in source}
        assert actual(target)
        assert quota_from_placement(graph, target) == source
    rng = Random(8811 + n)
    supports = [set(range(n)) - {v} for v in rng.sample(range(n), 3)]
    supports += [{v for v in range(n) if rng.random() < .8} for _ in range(6)]
    reverses = 0
    for selected in supports:
        feasible = actual(selected)
        if len(selected) == n - 1:
            assert feasible
        if feasible:
            recovered = quota_from_placement(graph, selected)
            assert len(recovered) <= len(selected) and source_valid(instance, recovered, local)
            round_truth(instance, recovered)
            reverses += 1
    print('Clause control:', m, 'clauses,', n, 'states,', len(lifts), 'lifts,',
          reverses, 'arbitrary reverse maps PASS', flush=True)


if __name__ == '__main__':
    local = clause_masks()
    exact_two_clause(local)
    for k in (1, 2, 4, 9):
        instance, selected = cycle_break(k)
        assert len(instance.clauses) == 8 * k + 2
        assert len(selected) == 7 * len(instance.clauses) + 1
        assert source_valid(instance, selected, local)
        for value in (False, True):
            assignment = {x: (value if x == 1 else True) for x in instance.ports}
            assert violated(instance.clauses, assignment) == k
        round_truth(instance, selected)
    print('Cyclic consistency loss at unbounded occurrence PASS', flush=True)
    sat = clause_quota_instance(((1, 2, 3), (-1, -2, -3)))
    unsat = clause_quota_instance(((1, 1, 1), (-1, -1, -1)))
    for instance in (sat, unsat):
        lifts = [lift_truth(instance, dict(zip(sorted(instance.ports), values)))
                 for values in product((False, True), repeat=len(instance.ports))]
        for length in (1, 3):
            check_graph(instance, lifts, local, length)
    instance, selected = cycle_break(1)
    check_graph(instance, [selected], local, 1)
