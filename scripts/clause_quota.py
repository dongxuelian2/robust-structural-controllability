"""Three-literal clauses in rank-three cotransversal matching blocks."""
from collections import defaultdict
from dataclasses import dataclass


CLAUSE_ROWS = (7, 25, 97, 84)


@dataclass
class ClauseQuota:
    clauses: tuple
    presentation: tuple
    elements: int
    groups: tuple
    forced: frozenset
    ports: dict


def clause_quota_instance(clauses):
    """Use a positive even number of signed three-literal clauses.

    Literal zero is invalid; repeated literals are permitted for tests.
    The standard hardness source has three distinct variables per clause.
    Elements 7c,...,7c+6 are p,a1,b1,a2,b2,a3,b3 of clause c.
    The last 3m elements are forced dummy resources.
    """
    clauses = tuple(tuple(clause) for clause in clauses)
    m = len(clauses)
    if not m or m % 2 or any(len(clause) != 3 or any(not x for x in clause)
                            for clause in clauses):
        raise ValueError('Use an even positive number of signed three-literal clauses.')
    rows = [row << (7 * c) for c in range(m) for row in CLAUSE_ROWS]
    for i in range(0, 3 * m, 2):
        rows.append(3 << (7 * m + i))
    ports = defaultdict(list)
    for c, clause in enumerate(clauses):
        for j, literal in enumerate(clause):
            a, b = 7 * c + 2 * j + 1, 7 * c + 2 * j + 2
            true, false = (a, b) if literal > 0 else (b, a)
            ports[abs(literal)].append((true, false, 7 * m + 3 * c + j))
    groups = [(7 * c, 7 * (c + 1)) for c in range(0, m, 2)]
    for occurrences in ports.values():
        for i, (true, _, dummy) in enumerate(occurrences):
            following_false = occurrences[(i + 1) % len(occurrences)][1]
            groups.append((true, following_false, dummy))
    forced = frozenset(range(0, 7 * m, 7)) | frozenset(range(7 * m, 10 * m))
    return ClauseQuota(clauses, tuple(rows), 10 * m, tuple(groups), forced, dict(ports))


def violated(clauses, assignment):
    return sum(not any(bool(assignment[abs(literal)]) == (literal > 0)
                       for literal in clause) for clause in clauses)


def lift_truth(instance, assignment):
    """A truth assignment violating u clauses gives cost 7m+u."""
    selected = set(instance.forced)
    for c, clause in enumerate(instance.clauses):
        satisfied = False
        for j, literal in enumerate(clause):
            truth = bool(assignment[abs(literal)]) == (literal > 0)
            selected.add(7 * c + 2 * j + (1 if truth else 2))
            satisfied |= truth
        if not satisfied:
            selected.add(7 * c + 1)
    assert len(selected) == 7 * len(instance.clauses) + violated(instance.clauses, assignment)
    return frozenset(selected)


def round_truth(instance, selected):
    """Round canonical occurrences by majority, within gamma(B)*e violations.

    Feasibility for matching and the group quotas is a precondition.
    The excess e is the cost above 7m. The loss gamma(B) is
    1+2*floor((B-1)/2), where B is maximum variable occurrence.
    """
    selected = set(selected)
    if not instance.forced <= selected or not selected <= set(range(instance.elements)):
        raise ValueError('A feasible placement must contain the forced resources.')
    for c in range(len(instance.clauses)):
        local = sum(7 * c + j in selected for j in range(7))
        if local < 4:
            raise ValueError('A clause matching block is not robust.')
    result = {}
    for variable, ports in instance.ports.items():
        values = [true in selected for true, false, _dummy in ports
                  if (true in selected) != (false in selected)]
        if len(values) == len(ports):
            assert len(set(values)) == 1
        result[variable] = 2 * sum(values) > len(values)
    excess = len(selected) - 7 * len(instance.clauses)
    bound = max(map(len, instance.ports.values()))
    loss = 1 + 2 * ((bound - 1) // 2)
    assert violated(instance.clauses, result) <= loss * excess
    return result
