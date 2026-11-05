"""Exact sparse control realization of uniform RAP with quota-two partitions."""
from collections import Counter
from dataclasses import dataclass
from assignment_transfer import robust_assignment, normalize_dominated
from bounded_core_solver import matching_number
from hardness_gadgets import protect_coloops
from scc_cover import pad_resource_network


@dataclass
class PartitionInstance:
    rows: tuple
    jobs: int
    options: tuple
    groups: tuple
    source_cycles: tuple
    base_rows: tuple
    neighborhoods: tuple
    projection: tuple


def partition_assignment_graph(n, options, groups, length=1):
    """Preserve OPT and all unit costs, with every state selectable.

    Original resources are partitioned into groups of size at least two,
    each requiring at least two selected resources. The original uniform
    robust assignment must be feasible. All neighborhoods are nonempty.
    The target is loopless, bipartite, has all undirected degrees in {2,3}
    and both directed degree bounds three, and has no mandatory location.
    Optional odd length gives the same optimum with multiplied girth.
    """
    options, groups = tuple(options), tuple(tuple(group) for group in groups)
    m = len(options)
    if (n < 1 or any(option <= 0 or option >> n for option in options) or
            any(len(group) < 2 for group in groups) or
            sorted(v for group in groups for v in group) != list(range(m)) or
            not robust_assignment(n, options, range(m))):
        raise ValueError('Use a feasible uniform RAP and a quota-two partition.')
    rows, neighborhoods = [0] * (2 * m + n), list(options)
    neighborhoods += [1 << j for j in range(n)] + [0] * m

    def attach(parent, members):
        if len(members) == 1:
            rows[parent] |= 1 << (m + members[0])
            return
        x, y = len(rows), len(rows) + 1
        rows.extend((1 << y, 1 << x))
        neighborhoods.extend((sum(1 << j for j in members), 0))
        rows[parent] |= 1 << x
        middle = len(members) // 2
        attach(y, members[:middle])
        attach(y, members[middle:])

    cycles = []
    for group in groups:
        cycle = []
        for i, e in enumerate(group):
            filler, following = m + n + e, group[(i + 1) % len(group)]
            rows[e] = 1 << filler
            rows[filler] = 1 << following
            attach(filler, [j for j in range(n) if options[following] >> j & 1])
            cycle.extend((e, filler))
        cycles.append(tuple(cycle))
    base = tuple(rows)
    projection = [v if option else -1 for v, option in enumerate(neighborhoods)]
    # Split only physical job columns, never original source resources.
    for v in range(m, m + n):
        while True:
            incoming = [u for u, row in enumerate(rows) if row >> v & 1]
            if len(incoming) <= 3:
                break
            x, y = len(rows), len(rows) + 1
            rows.extend((1 << y, (1 << x) | (1 << v)))
            neighborhoods.extend((neighborhoods[v], 0))
            projection.extend((v, -1))
            for u in incoming[:2]:
                rows[u] ^= (1 << v) | (1 << x)
    rows, cycles, neighborhoods, projection, _ = pad_resource_network(
        tuple(rows), tuple(cycles), tuple(neighborhoods), tuple(projection), {}, length)
    return PartitionInstance(rows, n, options, groups, cycles, base,
                             neighborhoods, projection)


def assignment_from_partition_placement(instance, selected):
    """Recover a no-larger robust assignment satisfying the original quotas.

    Feasibility of the input placement is a precondition. Rank projection,
    coloop protection and dominated-resource normalization are constructive;
    selected original resources are never removed at any step.
    """
    selected = set(selected)
    if not selected <= set(range(len(instance.rows))):
        raise ValueError('Invalid state index.')
    original_size = len(selected)
    for group, cycle in zip(instance.groups, instance.source_cycles):
        group, cycle = set(group), set(cycle)
        if len(selected & cycle) < 2:
            raise ValueError('A source quota is not satisfied.')
        assert all(instance.neighborhoods[v] == 0 for v in cycle - group)
        selected.difference_update(cycle - group)
        for v in sorted(group - selected):
            if len(selected & group) >= 2:
                break
            selected.add(v)
    assert len(selected) <= original_size
    base, n = instance.base_rows, instance.jobs
    full, nu = (1 << len(base)) - 1, len(base) - n
    active = tuple(v for v in range(len(base)) if instance.neighborhoods[v])
    index = {v: j for j, v in enumerate(active)}
    counts = Counter(instance.projection[v] for v in selected if instance.projection[v] in index)

    def rank(support):
        return len(support) + matching_number(
            base, full ^ sum(1 << v for v in support)) - nu

    normalized = protect_coloops(rank, active, counts)
    assert len(normalized) <= len(selected)
    result = normalize_dominated(n, instance.options,
                                tuple(instance.neighborhoods[v] for v in active),
                                {index[v] for v in normalized})
    assert len(result) <= original_size
    assert all(len(result & set(group)) >= 2 for group in instance.groups)
    return result
