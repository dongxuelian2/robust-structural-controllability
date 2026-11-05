"""Check the shared backup on published equation gadgets and control targets."""
from fractions import Fraction
from itertools import combinations, product
from random import Random

from assignment_transfer import assignment_from_placement, assignment_graph, robust_assignment
from bounded_core_solver import matching_number, source_components
from hardness_gadgets import complete_leaves, subdivide_ports
from shared_backup_cover import (
    cover_from_assignment, cover_from_quota_placement, shared_backup_instance, shared_backup_quota_graph,
)
from verify_quota_cotransversal import input_first_matchable
from verify_scc_cover import girth


def equation_gadget(parity, start=0):
    """Chlebík–Chlebíková Figure 2(i), including two contact-pair edges."""
    x0, x1, y0, y1, z0, z1, a0, a1 = (start + (i ^ parity) for i in range(8))
    s000, s011, s101, s110 = range(start + 8, start + 12)
    triangles = ((a0, s101, s110), (a1, s000, s011))
    packed_strong = ((x0, x1), (y0, s011), (y1, s000), (z0, s101), (z1, s110))
    packed_general = ((x1, a0), (x0, a1), *packed_strong[1:])
    edges = tuple(edge for triangle in triangles for edge in combinations(triangle, 2))
    edges += packed_general + ((x0, x1), (y0, y1), (a0, a1),
                               (y0, s110), (y1, s101), (z0, s011), (z1, s000))
    assert len(edges) == 17 + 2
    contacts = ((x0, x1), (y0, y1), (z0, z1))
    covers = []
    special = {(0, 0, 0): s000, (0, 1, 1): s011, (1, 0, 1): s101, (1, 1, 0): s110}
    for local in product((0, 1), repeat=3):
        independent = {pair[value] for pair, value in zip(contacts, local)}
        independent.add((a0, a1)[local[0]])
        if local in special:
            independent.add(special[local])
        covers.append(set(range(start, start + 12)) - independent)
    return edges, triangles, packed_general, packed_strong, contacts, covers


def covers_edges(edges, vertices):
    return all(u in vertices or v in vertices for u, v in edges)


def check_local():
    for parity in (0, 1):
        edges, triangles, general, strong, contacts, canonical = equation_gadget(parity)
        for is_strong, packed in ((False, general), (True, strong)):
            source = shared_backup_instance(12, edges, packed, triangles, is_strong)
            full = frozenset(range(len(source.options)))
            assert robust_assignment(source.jobs, source.options, full)
            for e in source.forced:
                selected = full - {e}
                assert (any(len(selected.intersection(g)) < 2 for g in source.groups)
                        or not robust_assignment(source.jobs, source.options, selected))
            count = 0
            for mask in range(1 << 12):
                cover = {v for v in range(12) if mask >> v & 1}
                selected = cover | source.forced
                feasible = (all(len(selected.intersection(g)) >= 2 for g in source.groups)
                            and robust_assignment(source.jobs, source.options, selected))
                is_cover = covers_edges(edges, cover)
                assert not feasible or is_cover
                if not is_strong or all(cover.intersection(pair) for pair in contacts):
                    assert feasible == is_cover
                count += feasible
            assert all(covers_edges(edges, cover) and
                       robust_assignment(source.jobs, source.options, cover | source.forced)
                       for cover in canonical)
            print('Equation', parity, 'strong', is_strong, '4096 supports,', count, 'feasible PASS', flush=True)


def check_control(strong):
    first, second = equation_gadget(0), equation_gadget(1, 12)
    edges = first[0] + second[0] + ((1, 12),)
    triangles = first[1] + second[1]
    index = 3 if strong else 2
    source = shared_backup_instance(24, edges, first[index] + second[index], triangles, strong)
    cover = first[-1][0] | second[-1][0]
    assert covers_edges(edges, cover)
    selected = cover | source.forced
    if strong:
        graph = subdivide_ports(complete_leaves(assignment_graph(source.jobs, source.options)))
        reps = range(len(graph.rows))
    else:
        graph = shared_backup_quota_graph(source, length=3)
        reps = graph.representatives
    rows = graph.rows
    n = len(rows)
    assert n - matching_number(rows, (1 << n) - 1) == source.jobs
    components = source_components(rows)
    if strong:
        assert components == [(1 << n) - 1]
    else:
        assert set(components) == {sum(1 << v for v in cycle) for cycle in graph.source_cycles}
    incoming = [sum(1 << u for u, row in enumerate(rows) if row >> v & 1) for v in range(n)]
    assert all(row and row.bit_count() <= 3 for row in incoming)
    assert all(not row >> v & 1 and 2 <= (row | incoming[v]).bit_count() <= 3
               for v, row in enumerate(rows))
    singleton_rows = [row for row in incoming if row.bit_count() == 1]
    assert len(singleton_rows) == len(set(singleton_rows))
    assert girth(rows) >= 12
    lifted = {reps[v] for v in selected}
    auxiliary = {reps[len(source.options)]} | set(Random(411 + strong).sample(range(n), 7))
    for placement in (lifted, lifted | auxiliary):
        mask = sum(1 << v for v in placement)
        assert all((mask & group).bit_count() >= 2 for group in components)
        assert all(input_first_matchable(rows, placement - {v}) for v in placement)
        if strong:
            recovered = assignment_from_placement(source.jobs, source.options, graph, placement)
            decoded = cover_from_assignment(source, recovered)
        else:
            decoded = cover_from_quota_placement(source, graph, placement)
        assert covers_edges(edges, decoded) and len(decoded) <= len(placement) - len(source.forced)
    print('Control strong', strong, n, 'states, girth', girth(rows), 'input failures and recovery PASS', flush=True)


if __name__ == '__main__':
    check_local()
    check_control(False)
    check_control(True)
    mu = Fraction(207, 10)
    assert 1 + Fraction(1, 2) / (4 * mu + 8) == Fraction(913, 908)
    assert 1 + Fraction(1, 2) / (7 * mu + 14) == Fraction(1594, 1589)
