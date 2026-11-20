"""Expander consistency in bicircular blocks with direct clause quotas."""
from collections import defaultdict
from dataclasses import dataclass
from math import isqrt
from assignment_transfer import branching_presentation, normalize_dominated


@dataclass(frozen=True)
class CopyBlock:
    start: int
    ports: int


@dataclass
class AmplifiedQuota:
    clauses: tuple
    predicate: str
    presentation: tuple
    elements: int
    groups: tuple
    forced: frozenset
    copies: dict
    degree: int
    options: dict
    relaxed: dict
    baseline: int


def copy_edges(ports):
    """Cycle, common edges, then two endpoint loops."""
    if ports not in (3, 4, 5, 6):
        raise ValueError('Use three through six paired ports.')
    order = 2 * ports
    true = [(v, v + 1) for v in range(0, order, 2)]
    false = [(v, (v + 1) % order) for v in range(1, order, 2)]
    common = [(2 * i + 1, 2 * i + 4) for i in range(ports - 2)]
    common.append((order - 3, order - 1))
    return tuple(true + false + common + [(0, 0), (2, 2)])


def copy_arcs(contacts, degree):
    """Eulerian directed expanders, specified by permutations.

    Degree four uses the four forward Margulis--Gabber--Galil maps.
    Degree two splits each vertex into a directed four-cycle with one
    incoming and outgoing external arc per port. A row lists its
    outgoing destinations in port order.
    """
    if contacts < 1 or degree not in (2, 4):
        raise ValueError('Use a positive contact count and degree two or four.')
    if contacts == 1:
        return ((0,) * degree,)
    cloud_size = 1 if degree == 4 else 4
    side = 1 + isqrt((contacts - 1) // cloud_size)
    base = []
    for x in range(side):
        for y in range(side):
            base.append((((x + 2 * y) % side) * side + y,
                         ((x + 2 * y + 1) % side) * side + y,
                         x * side + (y + 2 * x) % side,
                         x * side + (y + 2 * x + 1) % side))
    if degree == 4:
        return tuple(base)
    return tuple((4 * destination + i, 4 * u + (i + 1) % 4)
                 for u, row in enumerate(base) for i, destination in enumerate(row))


def amplified_quota_instance(clauses, degree=4, predicate='or'):
    """Build copy blocks and quota groups for OR or x OR (y XOR z)."""
    clauses = tuple(tuple(clause) for clause in clauses)
    m = len(clauses)
    if (degree not in (2, 4) or predicate not in ('or', 'oxr') or not m or m % 2 or
            any(len(clause) != 3 or any(not x for x in clause) for clause in clauses)):
        raise ValueError('Use an even positive clause count, OR or OXR, and degree two or four.')
    options, presentations, relaxed = {}, {}, {}
    for ports in range(degree + 1, degree + (3 if predicate == 'oxr' else 2)):
        options[ports] = tuple((1 << u) | (1 << v) for u, v in copy_edges(ports))
        presentations[ports], relaxed[ports] = branching_presentation(2 * ports, options[ports])
    occurrences = defaultdict(list)
    for c, clause in enumerate(clauses):
        for j, literal in enumerate(clause):
            occurrences[abs(literal)].append(3 * c + j)
    rows, groups, forced = [], [], set()
    copies, contacts = {}, [None] * (3 * m)
    size = baseline = 0
    for variable, positions in sorted(occurrences.items()):
        arcs = copy_arcs(len(positions), degree)
        blocks = []
        for i in range(len(arcs)):
            position = positions[i] if i < len(positions) else None
            special = predicate == 'oxr' and position is not None and position % 3 == 0
            ports = degree + 1 + int(special)
            block = CopyBlock(size, ports)
            blocks.append(block)
            if position is not None:
                contacts[position] = block
            rows.extend(row << size for row in presentations[ports])
            anchors = (size + 3 * ports - 1, size + 3 * ports)
            groups.append(anchors)
            forced.update(anchors)
            size += len(relaxed[ports])
            baseline += 2 * ports + 1
        copies[variable] = tuple(blocks)
        for u, outgoing in enumerate(arcs):
            for i, v in enumerate(outgoing):
                tail, head = blocks[u], blocks[v]
                groups.append((tail.start + i, head.start + head.ports + i,
                               tail.start + 2 * tail.ports + i))

    def literal_port(block, literal, index):
        return block.start + index + (0 if literal > 0 else block.ports)

    for c, clause in enumerate(clauses):
        blocks = contacts[3 * c:3 * c + 3]
        positive = tuple(literal_port(block, literal, degree) for block, literal in zip(blocks, clause))
        groups.append(positive + (size + c,))
        if predicate == 'oxr':
            x, y, z = clause
            second = (literal_port(blocks[0], x, degree + 1),
                      literal_port(blocks[1], -y, degree), literal_port(blocks[2], -z, degree))
            groups.append(second + (blocks[0].start + 2 * blocks[0].ports + degree,))
    rows.extend(3 << (size + i) for i in range(0, m, 2))
    forced.update(range(size, size + m))
    return AmplifiedQuota(clauses, predicate, tuple(rows), size + m, tuple(groups),
                          frozenset(forced), copies, degree, options, relaxed, baseline + m)


def copy_selection(ports, value):
    start = 0 if value else ports
    return set(range(start, start + ports)) | set(range(2 * ports, 3 * ports + 1))


def violations(instance, assignment):
    count = 0
    for clause in instance.clauses:
        x, y, z = [bool(assignment[abs(lit)]) == (lit > 0) for lit in clause]
        satisfied = x or (y != z) if instance.predicate == 'oxr' else x or y or z
        count += not satisfied
    return count


def lift_amplified(instance, assignment):
    """A truth assignment with u violations lifts at baseline plus u."""
    selected = set(instance.forced)
    for variable, blocks in instance.copies.items():
        for block in blocks:
            selected.update(block.start + e for e in copy_selection(block.ports, assignment[variable]))
    for group in instance.groups:
        if sum(e in selected for e in group) < 2:
            assert len(group) == 4 and group[-1] in selected
            selected.add(group[0])
    assert len(selected) == instance.baseline + violations(instance, assignment)
    return frozenset(selected)


def normalize_copy_resources(instance, selected):
    """Remove the added private job resources at no greater cost."""
    selected = set(selected)
    for blocks in instance.copies.values():
        for block in blocks:
            options, relaxed = instance.options[block.ports], instance.relaxed[block.ports]
            local = {e for e in range(len(relaxed)) if block.start + e in selected}
            if not any(e >= len(options) for e in local):
                continue
            recovered = normalize_dominated(2 * block.ports, options, relaxed, local)
            selected.difference_update(block.start + e for e in local - recovered)
            selected.update(block.start + e for e in recovered)
    return frozenset(selected)


def round_amplified(instance, selected):
    """Recover truth values from a feasible abstract placement."""
    normalized = normalize_copy_resources(instance, selected)
    assert instance.forced <= normalized and len(normalized) <= len(selected)
    assignment = {}
    for variable, blocks in instance.copies.items():
        votes = []
        for block in blocks:
            local = {e for e in range(len(instance.options[block.ports])) if block.start + e in normalized}
            minimum = 2 * block.ports + 1
            assert len(local) >= minimum
            if len(local) == minimum:
                value = 0 in local
                assert local == copy_selection(block.ports, value)
                votes.append(value)
        assignment[variable] = 2 * sum(votes) > len(votes)
    loss = 19 if instance.degree == 4 else 185
    assert violations(instance, assignment) <= loss * (len(normalized) - instance.baseline)
    return assignment
