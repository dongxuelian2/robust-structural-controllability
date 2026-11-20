"""Exact checks of optional-guard placements and the ternary row bound."""
from itertools import combinations, product

from assignment_transfer import assignment_graph, robust_assignment
from bounded_core_solver import matching_number, source_components
from direct_lin2 import parity_clause_jobs
from hardness_gadgets import complete_leaves, subdivide_ports
from shared_backup_cover import guarded_backup_assignment, shared_backup_quota_graph
from verify_direct_lin2 import build, EDGES
from verify_quota_cotransversal import input_first_matchable


def check_lift(values, clauses, buffers, optional, strong):
    choices = tuple((2 * v, 2 * v + 1) for v in range(len(values)))
    source = guarded_backup_assignment(2 * len(values), clauses, choices,
                                       strong=strong, optional=optional)
    selected = set(range(source.vertices, len(source.options)))
    selected.update(2 * v + value for v, value in enumerate(values))
    for v, j in buffers.items():
        assert j in optional
        selected.remove(source.options.index(1 << j, source.vertices + 2))
        selected.add(2 * v + 1 - values[v])
    assert len(selected) == source.jobs + 1
    assert source.forced <= selected
    assert all(len(selected.intersection(group)) >= 2 for group in source.groups)
    assert robust_assignment(source.jobs, source.options, selected)
    return source, selected


def check_plain(strong):
    triples = ((0, 1, 2), (0, 3, 4), (1, 3, 4))
    equations = tuple((*triple, rhs) for triple in triples for rhs in (0, 1))
    clauses = parity_clause_jobs(equations)
    optional = set(range(len(clauses)))
    for values in product((0, 1), repeat=5):
        false = {i: eq[:3] for i, eq in enumerate(equations)
                 if sum(values[v] for v in eq[:3]) % 2 != eq[3]}
        doubled = set(range(5))
        for v in range(5):
            if all((doubled - {v}).intersection(triple) for triple in false.values()):
                doubled.remove(v)
        buffers = {}
        for v in doubled:
            i, triple = next((i, triple) for i, triple in false.items()
                             if doubled.intersection(triple) == {v})
            pattern = tuple(2 * w + (1 - values[w] if w == v else values[w]) for w in triple)
            buffers[v] = clauses.index(pattern, 4 * i, 4 * i + 4)
        source, selected = check_lift(values, clauses, buffers, optional, strong)
    if strong:
        graph = subdivide_ports(complete_leaves(assignment_graph(source.jobs, source.options)))
        representatives = range(len(graph.rows))
    else:
        graph = shared_backup_quota_graph(source, length=3)
        representatives = graph.representatives
    placement = {representatives[v] for v in selected}
    states = len(graph.rows)
    mask = sum(1 << v for v in placement)
    assert all((mask & group).bit_count() >= 2 for group in source_components(graph.rows))
    assert states - matching_number(graph.rows, (1 << states) - 1) == source.jobs
    assert all(input_first_matchable(graph.rows, placement - {v}) for v in placement)
    print('Plain ternary strong', strong, '32 lifts and', states, 'state realization PASS', flush=True)


def check_amplified(strong):
    equations = ((0, 1, 2, 0), (0, 1, 2, 1), (0, 1, 2, 0))
    instance = build(equations, strong)
    stop = 2 * len(EDGES) * 3
    clauses = tuple(tuple(v for v in range(instance.source.vertices)
                          if instance.source.options[v] >> j & 1)
                    for j in range(stop + 4 * len(equations)))
    pairs = tuple(combinations(range(4), 2))
    for assignment in product((0, 1), repeat=3):
        values = tuple(assignment[v // instance.block_size] for v in range(len(instance.neighbors)))
        canonical = {2 * v + value for v, value in enumerate(values)}
        false = {i: eq[:3] for i, eq in enumerate(instance.equations)
                 if sum(values[v] for v in eq[:3]) % 2 != eq[3]}
        for start in range(len(pairs)):
            optional = {stop + 4 * i + j for i in range(len(equations))
                        for j in pairs[(start + i) % len(pairs)]}
            buffers = {}
            for i in false:
                job = next(j for j in range(stop + 4 * i, stop + 4 * (i + 1))
                           if j in optional and canonical.intersection(clauses[j]))
                spare, = set(clauses[job]) - canonical
                buffers[spare // 2] = job
            check_lift(values, clauses, buffers, optional, strong)
        optional = set(range(stop))
        buffers = {}
        for v, *_ in false.values():
            spare = 2 * v + 1 - values[v]
            buffers[v] = next(j for j in range(stop)
                              if spare in clauses[j] and canonical.intersection(clauses[j]))
        check_lift(values, clauses, buffers, optional, strong)
    print('Amplified ternary strong', strong, '56 lifts PASS', flush=True)


def check_rows(triples):
    variables = 3 * triples
    satisfying = [values for values in product((0, 1), repeat=variables)
                  if all(sum(values[i:i + 3]) % 2 == 0 for i in range(0, variables, 3))]
    largest = 0
    for signs in product((-1, 0, 1), repeat=variables):
        count = sum(sum(bit == sign for bit, sign in zip(values, signs)) == 1
                    for values in satisfying)
        largest = max(largest, count)
    assert 4 * largest == 3 * len(satisfying)
    print(triples, 'disjoint triples, sharp row count', largest, '/', len(satisfying), 'PASS')


if __name__ == '__main__':
    for strong in (False, True):
        check_plain(strong)
        check_amplified(strong)
    for triples in (1, 2):
        check_rows(triples)
