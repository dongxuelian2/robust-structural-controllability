"""Finite certificates for the expander clause construction."""
from itertools import product
from random import Random
from amplified_clause_quota import (
    amplified_quota_instance, copy_arcs, copy_edges, copy_selection,
    lift_amplified, normalize_copy_resources, round_amplified, violations,
)
from assignment_transfer import branching_presentation, normalize_dominated, robust_assignment
from bounded_core_solver import matching_number, source_components
from quota_cotransversal import quota_cotransversal_graph, quota_from_placement
from row_splitting import dual_rank_oracle
from verify_clause_quota import bitmask
from verify_gammoid_obstructions import robust, transversal_ranks
from verify_quota_cotransversal import input_first_matchable


def check_copy(ports):
    edges = copy_edges(ports)
    jobs = 2 * ports
    options = tuple((1 << u) | (1 << v) for u, v in edges)
    neighbors = [set() for _ in range(jobs)]
    for u, v in edges[:-2]:
        neighbors[u].add(v)
        neighbors[v].add(u)
    paths = []

    def walk(path, unused):
        if not unused:
            if 2 in neighbors[path[-1]]:
                paths.append(path + (2,))
            return
        for v in sorted(neighbors[path[-1]] & unused):
            walk(path + (v,), unused - {v})

    walk((0,), set(range(jobs)) - {0, 2})
    assert len(paths) == 2
    rows = [sum(1 << e for e, mask in enumerate(options) if mask >> j & 1) for j in range(jobs)]
    ranks = transversal_ranks(rows, len(options))
    anchored = 3 << (len(options) - 2)
    minimum = {mask for mask in range(1 << len(options))
               if mask & anchored == anchored and mask.bit_count() == jobs + 1 and robust(ranks, mask)}
    assert minimum == {bitmask(copy_selection(ports, value)) for value in (False, True)}
    presentation, relaxed = branching_presentation(jobs, options)
    dual = dual_rank_oracle(presentation, len(relaxed))
    rng = Random(1979 + ports)
    exchanges = 0
    for _ in range(192):
        selected = {e for e in range(len(relaxed)) if rng.random() < .75}
        originals = selected & set(range(len(options)))
        assert dual(originals) == ranks[bitmask(originals)]
        if not robust_assignment(jobs, relaxed, selected):
            continue
        recovered = normalize_dominated(jobs, options, relaxed, selected)
        assert originals <= recovered and len(recovered) <= len(selected)
        assert robust(ranks, bitmask(recovered))
        exchanges += 1
    print('Copy block:', jobs, 'jobs, paths', paths, ';', exchanges, 'dominated exchanges PASS', flush=True)


def check_arcs(degree):
    audited = set()
    cuts = 0
    for contacts in range(1, 65):
        arcs = copy_arcs(contacts, degree)
        n = len(arcs)
        assert contacts <= n <= (2 if degree == 4 else 4) * contacts
        assert all(sorted(row[i] for row in arcs) == list(range(n)) for i in range(degree))
        if n in audited or n > 16:
            continue
        audited.add(n)
        for mask in range(1, 1 << n):
            if mask.bit_count() > n // 2:
                continue
            boundary = sum(not (mask >> v & 1) for u, row in enumerate(arcs) if mask >> u & 1 for v in row)
            if degree == 4:
                assert 9 * boundary >= 2 * mask.bit_count()
            else:
                assert 92 * boundary >= mask.bit_count()
            cuts += 1
    print('Directed consistency degree:', degree, ';', cuts, 'exact small cuts PASS', flush=True)


def abstract_feasible(instance, selected):
    if not instance.forced <= selected or any(len(selected & set(group)) < 2 for group in instance.groups):
        return False
    for copies in instance.copies.values():
        for block in copies:
            relaxed = instance.relaxed[block.ports]
            local = {e for e in range(len(relaxed)) if block.start + e in selected}
            if not robust_assignment(2 * block.ports, relaxed, local):
                return False
    return True


def mixed_support(instance):
    """Make a feasible support that uses different labels and repairs copies."""
    selected = set(instance.forced)
    owner = {block.start + e: block for copies in instance.copies.values()
             for block in copies for e in range(len(instance.options[block.ports]))}
    for copies in instance.copies.values():
        for i, block in enumerate(copies):
            selected.update(block.start + e for e in copy_selection(block.ports, i % 2 == 0))
    repaired = set()
    for group in instance.groups:
        if len(selected & set(group)) >= 2:
            continue
        # For an implication use its tail, for a clause use its first literal.
        repaired.add(owner[group[0]])
    for block in repaired:
        count = len(instance.options[block.ports])
        selected.update(range(block.start, block.start + count + 1))
    return frozenset(selected), len(repaired)


def check_formula(degree, predicate):
    fixtures = [((1, 2, 3), (-1, -2, -3)), ((1, 1, 1), (-1, -1, -1)),
                ((1, 2, 3), (-1, 2, 4), (1, -2, 3), (-1, -2, -4))]
    forward = reverse = 0
    for formula in fixtures:
        instance = amplified_quota_instance(formula, degree, predicate)
        grouped = [e for group in instance.groups for e in group]
        assert len(grouped) == len(set(grouped))
        for values in product((False, True), repeat=len(instance.copies)):
            truth = dict(zip(instance.copies, values))
            source = lift_amplified(instance, truth)
            assert abstract_feasible(instance, source)
            recovered = round_amplified(instance, source)
            if not violations(instance, truth):
                assert recovered == truth
            forward += 1
        selected, repaired = mixed_support(instance)
        assert repaired and abstract_feasible(instance, selected)
        normalized = normalize_copy_resources(instance, selected)
        assert len(normalized) < len(selected) and abstract_feasible(instance, normalized)
        round_amplified(instance, selected)
        reverse += 1
    print('Direct clause quotas:', predicate, degree, ';', forward, 'truth lifts,',
          reverse, 'repaired supports with private resources PASS', flush=True)


def check_control(degree, length, predicate):
    instance = amplified_quota_instance(((1, 2, 3), (-4, 5, 6)), degree, predicate)
    target = quota_cotransversal_graph(instance.presentation, instance.elements, instance.groups, length)
    rows, size = target.rows, len(target.rows)
    full = (1 << size) - 1
    dual = dual_rank_oracle(instance.presentation, instance.elements)
    assert size - matching_number(rows, full) == dual(range(instance.elements))
    sources = source_components(rows)
    assert set(sources) == {bitmask(cycle) for cycle in target.source_cycles}
    source = lift_amplified(instance, {x: x == 1 for x in instance.copies})
    assert len(source) == instance.baseline
    lifted = {target.representatives[e] for e in source}
    block = next(iter(instance.copies.values()))[0]
    added = target.representatives[block.start + len(instance.options[block.ports])]
    # Choose a private resource and actual auxiliary states in a feasible support.
    extras = set(Random(551 + degree).sample(range(size), 5)) | {added}
    for selected in (lifted, lifted | extras):
        assert all((bitmask(selected) & group).bit_count() >= 2 for group in sources)
        assert all(input_first_matchable(rows, selected - {e}) for e in selected)
        recovered = quota_from_placement(target, selected)
        assert len(recovered) <= len(selected) and abstract_feasible(instance, recovered)
        assert violations(instance, round_amplified(instance, recovered)) == 0
    print('Amplified control:', predicate, degree, '; length', length, ';', size,
          'states, direct input failures and recovery PASS', flush=True)


if __name__ == '__main__':
    for ports in (3, 4, 5, 6):
        check_copy(ports)
    for degree in (2, 4):
        check_arcs(degree)
        for predicate in ('or', 'oxr'):
            check_formula(degree, predicate)
    check_control(4, 1, 'oxr')
    check_control(2, 3, 'or')
