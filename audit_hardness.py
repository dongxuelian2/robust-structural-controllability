"""Fresh audits of the subdivision lemma and its constructive reverse map.

The universal small-presentation test enumerates matchings as injections,
without the existing rank or matching routines. The normalization audit
uses explicit dedicated-input rows and tests every failure of the support.
Python 3.10+, standard library only.
"""
from itertools import combinations, permutations, product
from hardness_gadgets import (vertex_cover_graph, subdivide_ports,
                              normalize_subdivision, extract_cover)


def transversal_bases(rows, right_count):
    return {sum(1 << v for v in assignment)
            for assignment in permutations(range(right_count), len(rows))
            if all(row >> v & 1 for row, v in zip(rows, assignment))}


def audit_edge_subdivision():
    presentations = operations = splits = 0
    for left_count in range(1, 4):
        for right_count in range(left_count, 5):
            full = (1 << right_count) - 1
            for rows in product(range(full + 1), repeat=left_count):
                bases = transversal_bases(rows, right_count)
                if not bases:
                    continue
                presentations += 1
                dual = {full ^ b for b in bases}
                for u, row in enumerate(rows):
                    for v in range(right_count):
                        if not row >> v & 1:
                            continue
                        alpha = 1 << right_count
                        modified = list(rows)
                        modified[u] = (row ^ (1 << v)) | alpha
                        modified.append((1 << v) | alpha)
                        new_full = full | alpha
                        actual = {new_full ^ b for b in
                                  transversal_bases(modified, right_count + 1)}
                        expected = dual | {(b ^ (1 << v)) | alpha
                                           for b in dual if b >> v & 1}
                        assert actual == expected, (rows, u, v, actual, expected)
                        operations += 1
                # Moving several incidences at once is the target-splitting
                # operation needed to reduce a high target degree.
                for v in range(right_count):
                    incident = [u for u, row in enumerate(rows) if row >> v & 1]
                    for subset in range(1, 1 << len(incident)):
                        if subset.bit_count() < 2:
                            continue  # already covered by the single-edge test
                        alpha = 1 << right_count
                        modified = list(rows)
                        for j, u in enumerate(incident):
                            if subset >> j & 1:
                                modified[u] = (modified[u] ^ (1 << v)) | alpha
                        modified.append((1 << v) | alpha)
                        actual = {(full | alpha) ^ b for b in
                                  transversal_bases(modified, right_count + 1)}
                        expected = dual | {(b ^ (1 << v)) | alpha
                                           for b in dual if b >> v & 1}
                        assert actual == expected, (rows, v, subset)
                        splits += 1
    print('Full-row-rank presentations:', presentations,
          '; individual edge subdivisions:', operations,
          '; multi-incidence target splits:', splits, 'PASS', flush=True)


def augmented_matchable(rows, selected):
    """Direct state-plus-dedicated-input matching; no dual-rank formula."""
    n = len(rows)
    all_rows = list(rows) + [1 << v for v in selected]
    mate = [-1] * n

    def augment(u, seen):
        remaining = all_rows[u]
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            v = bit.bit_length() - 1
            if v in seen:
                continue
            seen.add(v)
            if mate[v] < 0 or augment(mate[v], seen):
                mate[v] = u
                return True
        return False

    return sum(augment(u, set()) for u in range(len(all_rows))) == n


def audit_normalization(edges):
    base = vertex_cover_graph(edges)
    enlarged = subdivide_ports(base)
    relevant = sorted(set(base.owner) |
                      {v for pair in base.forced_pairs for v in pair})
    copies = {v: [] for v in relevant}
    for v, origin in enumerate(enlarged.projection):
        if origin in copies:
            copies[origin].append(v)
    assert all(len(c) >= 2 for c in copies.values())
    count = len(relevant)
    ordinary = [augmented_matchable(base.rows,
                [v for j, v in enumerate(relevant) if mask >> j & 1])
                for mask in range(1 << count)]
    robust_sets = [ordinary[mask] and all(ordinary[mask ^ (1 << j)]
                   for j in range(count) if mask >> j & 1)
                   for mask in range(1 << count)]
    index = {v: j for j, v in enumerate(relevant)}
    checked = feasible = 0
    for multiplicity in product(range(3), repeat=count):
        checked += 1
        support = sum(1 << j for j, amount in enumerate(multiplicity) if amount)
        if not ordinary[support] or any(not ordinary[support ^ (1 << j)]
                for j, amount in enumerate(multiplicity) if amount == 1):
            continue
        feasible += 1
        selected = {copy for v, amount in zip(relevant, multiplicity)
                    for copy in copies[v][:amount]}
        normalized = normalize_subdivision(base, enlarged, selected)
        mask = sum(1 << index[v] for v in normalized)
        assert robust_sets[mask], (edges, multiplicity, normalized)
        assert len(normalized) <= sum(multiplicity)
        decoded = extract_cover(base, normalized)
        assert all(set(e) & decoded for e in edges)
        assert len(decoded) <= len(selected) - 2 * len(edges)
    print('Normalization source', edges, ':', checked, 'multiplicities,',
          feasible, 'robust placements PASS', flush=True)


def euler_characteristic(rotation):
    """Count faces in an orientable rotation system, using only integer darts."""
    darts = {(u, v) for u, neighbors in rotation.items() for v in neighbors}
    assert all((v, u) in darts for u, v in darts)
    unused, faces = set(darts), 0
    while unused:
        start = next(iter(unused))
        current = start
        while True:
            unused.remove(current)
            u, v = current
            neighbors = rotation[v]
            current = v, neighbors[(neighbors.index(u) + 1) % len(neighbors)]
            if current == start:
                break
        faces += 1
    return len(rotation) - len(darts) // 2 + faces


def port_rotation(instance):
    """Expand a source rotation into the port construction's rotation system."""
    index = {label: v for v, label in enumerate(instance.labels)}
    result = {u: [v for v in range(len(instance.rows)) if row >> v & 1]
              for u, row in enumerate(instance.rows)}
    for u, label in enumerate(instance.labels):
        if label[0] != 'p':
            continue
        _, vertex, position = label
        forward = index['q', vertex, position]
        previous = index.get(('q', vertex, position - 1))
        if previous is None:
            previous = index.get(('prefix', vertex))
        element = next(v for v in result[u] if instance.labels[v][0] == 'x')
        result[u] = ([] if previous is None else [previous]) + [element, forward]
    assert all(set(result[u]) == {v for v in range(len(instance.rows))
                                 if instance.rows[u] >> v & 1} for u in result)
    return result


def audit_planarity_ports():
    # A planar wheel whose default lexicographic ports are not its rotation.
    edges = [(0, v) for v in range(1, 5)] + [(1, 3), (3, 2), (2, 4), (4, 1)]
    source_rotation = {0: [1, 4, 2, 3], 1: [0, 3, 4], 2: [3, 0, 4],
                       3: [1, 0, 2], 4: [2, 0, 1]}
    assert euler_characteristic(source_rotation) == 2
    wrong = vertex_cover_graph(edges)
    right = vertex_cover_graph(edges, port_order=source_rotation)
    # Nine internally disjoint paths certify K_3,3 in the default target.
    paths = [[34, 7, 46, 47, 48], [34, 33, 32, 3, 40, 41, 42], [34, 35, 36],
             [36, 11, 52, 53, 54], [36, 37, 38, 15, 58, 59, 60],
             [42, 19, 54], [42, 43, 44, 23, 60], [48, 27, 56, 55, 54],
             [48, 49, 50, 31, 62, 61, 60]]
    sides = ({34, 54, 60}, {36, 42, 48})
    branches = sides[0] | sides[1]
    seen_internal, pairs = set(), set()
    for path in paths:
        assert all(wrong.rows[u] >> v & 1 for u, v in zip(path, path[1:]))
        assert not set(path[1:-1]) & (seen_internal | branches)
        seen_internal.update(path[1:-1])
        u, v = path[0], path[-1]
        if u in sides[1]:
            u, v = v, u
        pairs.add((u, v))
    assert pairs == {(u, v) for u in sides[0] for v in sides[1]}
    assert euler_characteristic(port_rotation(wrong)) != 2
    assert euler_characteristic(port_rotation(right)) == 2
    print('Wheel: default K_3,3 certificate; corrected genus-zero rotation PASS')


def audit_raw_incidence():
    """Bicircular duals: a complementary spanning tree is always optimal."""
    count = 0
    for n in (4, 5):
        pairs = list(combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            edges = [e for j, e in enumerate(pairs) if mask >> j & 1]
            if len(edges) <= n:
                continue
            seen, tree = {0}, set()
            while True:
                crossing = next((e for e in edges if
                                 (e[0] in seen) != (e[1] in seen)), None)
                if crossing is None:
                    break
                seen.update(crossing)
                tree.add(crossing)
            if len(seen) != n:
                continue
            rows = [0] * (n + len(edges))
            for j, (u, v) in enumerate(edges):
                for w in (u, v):
                    rows[w] |= 1 << (n + j)
                    rows[n + j] |= 1 << w
            selected = {n + j for j, e in enumerate(edges) if e not in tree}
            assert len(selected) == len(edges) - n + 1
            assert augmented_matchable(rows, selected)
            assert all(augmented_matchable(rows, selected - {v}) for v in selected)
            # At most 2n state edges can be matched in this bipartite graph.
            # A failure leaves m-n inputs, so these full matchings also
            # certify state matching size 2n and the lower bound m-n+1.
            count += 1
    print('Raw bicircular incidence collapse:', count, 'connected sources PASS')


if __name__ == '__main__':
    audit_edge_subdivision()
    for edges in (((0, 1),), ((0, 1), (1, 2)), ((0, 1), (2, 3))):
        audit_normalization(edges)
    audit_planarity_ports()
    audit_raw_incidence()
