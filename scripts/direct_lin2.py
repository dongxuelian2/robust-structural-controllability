"""Shared-backup reductions from three-variable parity systems."""
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import product

from shared_backup_cover import guarded_backup_assignment


def parity_clause_jobs(equations):
    """Return four literal neighborhoods for each equation."""
    return tuple(tuple(2 * x + 1 - sign for x, sign in zip((u, v, w), signs))
                 for u, v, w, rhs in equations for signs in product((0, 1), repeat=3)
                 if sum(signs) % 2 != rhs)


def parity_deletion_instance(variables, equations, strong=False):
    """Transfer the minimum number of variables to delete for consistency."""
    equations = tuple(equations)
    if (variables < 1 or any(len(eq) != 4 or eq[3] not in (0, 1)
                            or len(set(eq[:3])) != 3
                            or any(not 0 <= x < variables for x in eq[:3])
                            for eq in equations)):
        raise ValueError('Use three distinct variables per equation.')
    packed = tuple((2 * v, 2 * v + 1) for v in range(variables))
    return guarded_backup_assignment(2 * variables, parity_clause_jobs(equations),
                                     packed, strong=strong)


@dataclass(frozen=True)
class Lin2Instance:
    block_size: int
    neighbors: tuple
    contacts: tuple
    equations: tuple
    source: object


def lin2_instance(block_size, edges, contacts, equations, strong=False):
    """Use the same amplifier for each variable.

    Variables are indexed by consecutive integers starting at zero.
    Each variable has one occurrence per amplifier contact.
    The amplifier cut property is assumed.
    """
    contacts, equations = tuple(contacts), tuple(equations)
    edges = tuple(tuple(sorted(edge)) for edge in edges)
    if (not contacts or len(set(contacts)) != len(contacts)
            or any(not 0 <= v < block_size for v in contacts)
            or len(set(edges)) != len(edges)
            or any(len(edge) != 2 or not 0 <= edge[0] < edge[1] < block_size for edge in edges)
            or any(len(eq) != 4 or eq[3] not in (0, 1) for eq in equations)):
        raise ValueError('Use a simple amplifier and three-variable parity equations.')
    degree = Counter(v for edge in edges for v in edge)
    if any(degree[v] != (2 if v in contacts else 3) for v in range(block_size)):
        raise ValueError('Contacts have degree two and checkers have degree three.')
    occurrences = defaultdict(list)
    local = [list(eq) for eq in equations]
    for i, equation in enumerate(equations):
        for j, x in enumerate(equation[:3]):
            occurrences[x].append((i, j))
    if (not occurrences or sorted(occurrences) != list(range(len(occurrences)))
            or any(len(positions) != len(contacts) for positions in occurrences.values())):
        raise ValueError('Every variable needs exactly one occurrence per contact.')
    size = block_size * len(occurrences)
    neighbors = [[] for _ in range(size)]
    private_jobs, contact_blocks = [], []
    for x, positions in sorted(occurrences.items()):
        start = block_size * x
        copies = tuple(start + v for v in contacts)
        contact_blocks.append(copies)
        for (i, j), v in zip(positions, copies):
            local[i][j] = v
        for u, v in edges:
            u, v = start + u, start + v
            neighbors[u].append(v)
            neighbors[v].append(u)
            private_jobs.extend(((2 * u, 2 * v + 1), (2 * u + 1, 2 * v)))
    private_jobs.extend(parity_clause_jobs(local))
    packed_jobs = tuple((2 * v, 2 * v + 1) for v in range(size))
    source = guarded_backup_assignment(2 * size, tuple(private_jobs), packed_jobs, strong=strong)
    return Lin2Instance(block_size, tuple(tuple(row) for row in neighbors),
                        tuple(contact_blocks), tuple(tuple(eq) for eq in local), source)


def violated(values, equations):
    return sum(values[u] ^ values[v] ^ values[w] != rhs for u, v, w, rhs in equations)


def lift_lin2(instance, assignment):
    values = [assignment[v // instance.block_size] for v in range(len(instance.neighbors))]
    selected = set(instance.source.forced)
    selected.update(2 * v + value for v, value in enumerate(values))
    for u, v, w, rhs in instance.equations:
        if values[u] ^ values[v] ^ values[w] != rhs:
            selected.add(2 * u + 1 - values[u])
    return frozenset(selected)


def round_lin2(instance, selected):
    """Round the selected original labels of a feasible placement."""
    size = len(instance.neighbors)
    values = [int(2 * v + 1 in selected) for v in range(size)]
    assert all(2 * v + value in selected for v, value in enumerate(values))
    doubled = [v for v in range(size) if 2 * v in selected and 2 * v + 1 in selected]
    equation_at = {v: eq for eq in instance.equations for v in eq[:3]}
    changed = True
    while changed:
        changed = False
        for v in doubled:
            neighbors = instance.neighbors[v]
            difference = len(neighbors) - 2 * sum(values[v] != values[u] for u in neighbors)
            if v in equation_at:
                difference += 1 - 2 * violated(values, (equation_at[v],))
            if difference < 0:
                values[v] ^= 1
                changed = True
    cut = sum(values[v] != values[u] for v, row in enumerate(instance.neighbors) for u in row) // 2
    local_loss = cut + violated(values, instance.equations)
    assert local_loss <= len(doubled)
    return tuple(int(2 * sum(values[v] for v in copies) > len(copies))
                 for copies in instance.contacts)
