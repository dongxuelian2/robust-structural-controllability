"""Check parity presentations and the exact variable-deletion transfer."""
from fractions import Fraction
from itertools import combinations, product

from assignment_transfer import assignment_from_placement, assignment_graph, robust_assignment
from bounded_core_solver import source_components
from direct_lin2 import parity_deletion_instance, violated
from hardness_gadgets import complete_leaves, subdivide_ports
from row_splitting import dual_rank_oracle
from shared_backup_cover import cover_from_quota_placement, shared_backup_quota_graph
from verify_quota_cotransversal import input_first_matchable


def check_columns(columns, label_count, expected):
    size = len(columns)
    rows = tuple(sum(1 << e for e, mask in enumerate(columns) if mask >> j & 1)
                 for j in range(max(columns).bit_length()))
    rank = dual_rank_oracle(rows, size)
    full_rank = rank(range(size))
    for omitted in range(1 << label_count):
        selected = {e for e in range(size) if not omitted >> e & 1}
        feasible = all(rank(selected - {e}) == full_rank for e in selected)
        assert feasible == expected(omitted)


def check_tight(arity):
    pairs = tuple(3 << (2 * i) for i in range(arity))
    labels = tuple(1 << v for v in range(2 * arity))
    bad = {sum(1 << (2 * i + bit) for i, bit in enumerate(bits))
           for bits in product((0, 1), repeat=arity) if sum(bits) % 2}

    def expected(mask):
        return mask not in bad and all(mask & pair != pair for pair in pairs)

    check_columns(labels + pairs + tuple(sorted(bad)), 2 * arity, expected)


def check_k4():
    # Opposite edges form pairs (01,23), (02,13), (03,12).
    columns = (3, 12, 5, 10, 9, 6, 1, 2)

    def expected(mask):
        count = mask.bit_count()
        if count < 3:
            return True
        if count > 3:
            return False
        # Three edges are either a spanning tree or a triangle.
        covered = 0
        for v in range(6):
            if mask >> v & 1:
                covered |= columns[v]
        return covered == 15

    check_columns(columns, 6, expected)
    canonical = [sum(1 << (2 * i + bit) for i, bit in enumerate(bits))
                 for bits in product((0, 1), repeat=3)]
    assert [expected(mask) for mask in canonical] == [True, False, False, True, False, True, True, False]
    assert expected(0b111)


def check_fill():
    equations = tuple((*edge, rhs) for edge in combinations(range(4), 3) for rhs in (0, 1))
    for partial in product((-1, 0, 1), repeat=4):
        deleted = [i for i, bit in enumerate(partial) if bit == -1]
        untouched = [eq for eq in equations if not set(eq[:3]).intersection(deleted)]
        full = []
        for bits in product((0, 1), repeat=len(deleted)):
            assignment = list(partial)
            for i, bit in zip(deleted, bits):
                assignment[i] = bit
            full.append(violated(assignment, equations))
        expected = violated(partial, untouched) + Fraction(len(equations) - len(untouched), 2)
        assert Fraction(sum(full), len(full)) == expected
    assert Fraction(7**7, 3**6) * (Fraction(5, 7) + Fraction(1, 10000))**21 < 1


def check_assignment(strong):
    equations = ((0, 1, 2, 0), (0, 1, 2, 1))
    source = parity_deletion_instance(3, equations, strong)
    minimum = len(source.options)
    for mask in range(1 << 6):
        labels = {v for v in range(6) if mask >> v & 1}
        selected = labels | source.forced
        pairs = [{2 * i, 2 * i + 1} & labels for i in range(3)]
        expected = all(pairs) and any(len(pair) == 2 for pair in pairs)
        assert robust_assignment(source.jobs, source.options, selected) == expected
        if expected:
            minimum = min(minimum, len(selected))
    baseline = 3 + (8 if strong else 4) * len(equations) + 2
    assert minimum == baseline + 1
    if strong:
        graph = subdivide_ports(complete_leaves(assignment_graph(source.jobs, source.options)))
        reps = range(len(graph.rows))
    else:
        graph = shared_backup_quota_graph(source, length=3)
        reps = graph.representatives
    lifted = {reps[v] for v in source.forced | {0, 1, 2, 4}}
    for selected in (lifted, lifted | {reps[len(source.options)]}):
        mask = sum(1 << v for v in selected)
        assert all((mask & group).bit_count() >= 2 for group in source_components(graph.rows))
        assert all(input_first_matchable(graph.rows, selected - {v}) for v in selected)
        if strong:
            original = assignment_from_placement(source.jobs, source.options, graph, selected)
            labels = {v for v in original if v < 6}
        else:
            labels = cover_from_quota_placement(source, graph, selected)
        assert len(labels) + len(source.forced) <= len(selected)
        assert all(labels.intersection((2 * i, 2 * i + 1)) for i in range(3))
        assert any({2 * i, 2 * i + 1} <= labels for i in range(3))
    print('Deletion strong', strong, '64 supports and', len(graph.rows), 'state control PASS')


if __name__ == '__main__':
    check_tight(3)
    check_tight(4)
    check_k4()
    check_fill()
    print('Parity presentations, 81 random-fill identities and exponent PASS')
    for strong in (False, True):
        check_assignment(strong)
