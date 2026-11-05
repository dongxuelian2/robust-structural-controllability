"""Exact checks of direct parity costs, rounding and actuator failures."""
from fractions import Fraction
from itertools import product
from random import Random

from assignment_transfer import assignment_from_placement, assignment_graph, robust_assignment
from bounded_core_solver import source_components
from direct_lin2 import lift_lin2, lin2_instance, round_lin2, violated
from hardness_gadgets import complete_leaves, subdivide_ports
from shared_backup_cover import cover_from_quota_placement, shared_backup_quota_graph
from verify_quota_cotransversal import input_first_matchable


CONTACTS = (2, 3, 4)
EDGES = tuple(product((0, 1), CONTACTS))


def build(equations, strong):
    return lin2_instance(5, EDGES, CONTACTS, equations, strong)


def check_amplifier():
    for mask in range(1 << 5):
        contacts = sum(mask >> v & 1 for v in CONTACTS)
        cut = sum((mask >> u ^ mask >> v) & 1 for u, v in EDGES)
        assert cut >= min(contacts, len(CONTACTS) - contacts)


def check_rounding(strong):
    for rhs in (0, 1):
        equations = ((0, 0, 0, rhs),)
        instance = build(equations, strong)
        source = instance.source
        # Include empty pairs and compare robustness with the job conditions.
        for mask in range(1 << source.vertices):
            labels = {v for v in range(source.vertices) if mask >> v & 1}
            selected = labels | source.forced
            all_choices = (1 << (2 * len(EDGES) + 4 + len(instance.neighbors))) - 1
            supported = 0
            for v in labels:
                supported |= source.options[v]
            feasible = robust_assignment(source.jobs, source.options, selected)
            assert feasible == (supported == all_choices)
            if feasible:
                loss = violated(round_lin2(instance, labels), equations)
                assert loss <= len(labels) - len(instance.neighbors)
    equations = ((0, 1, 2, 0), (0, 1, 2, 1), (0, 1, 2, 0))
    instance = build(equations, strong)
    source = instance.source
    baseline = len(source.forced) + len(instance.neighbors)
    for values in product((0, 1), repeat=3):
        selected = lift_lin2(instance, values)
        assert len(selected) == baseline + violated(values, equations)
        assert robust_assignment(source.jobs, source.options, selected)
    rng = Random(571)
    for _ in range(100):
        # Choose mixed labels and repair each deficient literal group.
        labels = {2 * v + rng.randrange(2) for v in range(len(instance.neighbors))}
        for job in range(2 * len(EDGES) * 3 + 4 * len(equations)):
            eligible = {v for v in range(source.vertices) if source.options[v] >> job & 1}
            if not labels.intersection(eligible):
                labels.add(min(eligible))
        selected = labels | source.forced
        assert robust_assignment(source.jobs, source.options, selected)
        assert 1 <= violated(round_lin2(instance, labels), equations) <= len(selected) - baseline
    print('Parity strong', strong, '2048 exact supports, 108 lifts and mixed supports PASS', flush=True)


def check_control(strong):
    equations = ((0, 0, 0, 1),)
    instance = build(equations, strong)
    source = instance.source
    if strong:
        graph = subdivide_ports(complete_leaves(assignment_graph(source.jobs, source.options)))
        reps = range(len(graph.rows))
    else:
        graph = shared_backup_quota_graph(source, length=3)
        reps = graph.representatives
    lifted = {reps[v] for v in lift_lin2(instance, (1,))}
    auxiliary = {reps[len(source.options)]} | set(Random(369).sample(range(len(graph.rows)), 5))
    for selected in (lifted, lifted | auxiliary):
        mask = sum(1 << v for v in selected)
        assert all((mask & group).bit_count() >= 2 for group in source_components(graph.rows))
        assert all(input_first_matchable(graph.rows, selected - {v}) for v in selected)
        if strong:
            original = assignment_from_placement(source.jobs, source.options, graph, selected)
            labels = {v for v in original if v < source.vertices}
        else:
            labels = cover_from_quota_placement(source, graph, selected)
        assert len(labels) <= len(selected) - len(source.forced)
        assert violated(round_lin2(instance, labels), equations) == 0
    print('Parity control strong', strong, len(graph.rows), 'states, input failures and recovery PASS', flush=True)


if __name__ == '__main__':
    check_amplifier()
    for strong in (False, True):
        check_rounding(strong)
        check_control(strong)
    mu = Fraction(207, 10)
    assert 1 + Fraction(1, 2) / (4 * mu + 1) == Fraction(843, 838)
    assert 1 + Fraction(1, 2) / (7 * mu + 2) == Fraction(1474, 1469)
