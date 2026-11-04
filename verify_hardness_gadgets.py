"""Independent integer matching checks for the hardness and girth gadgets."""
from itertools import combinations
from random import Random
from hardness_gadgets import (vertex_cover_graph, lift_cover, extract_cover,
                              subdivide_ports, complete_leaves, normalize_subdivision)
from verify_outdegree_two import matching_subsets, brute_robust, independent_sources


def match(rows, allowed):
    """Recursive augmenting oracle, independent of the repository's BFS matcher."""
    n = len(rows)
    mate = [-1] * n

    def augment(u, seen):
        available = rows[u] & allowed
        while available:
            bit = available & -available
            available ^= bit
            v = bit.bit_length() - 1
            if v in seen:
                continue
            seen.add(v)
            if mate[v] < 0 or augment(mate[v], seen):
                mate[v] = u
                return True
        return False

    return sum(augment(v, set()) for v in range(n))


def sc(rows, selected):
    full = (1 << len(rows)) - 1
    mask = sum(1 << v for v in selected)
    return bool(mask) and match(rows, full ^ mask) == len(rows) - len(selected)


def robust(rows, selected):
    selected = set(selected)
    return sc(rows, selected) and all(sc(rows, selected - {v}) for v in selected)


def source_covers(edges):
    vertices = sorted({v for e in edges for v in e})
    covers = [set(c) for size in range(len(vertices) + 1)
              for c in combinations(vertices, size)
              if all(set(e) & set(c) for e in edges)]
    return covers


def structural_checks(instance):
    rows = instance.rows
    n = len(rows)
    assert all(not row >> v & 1 and row.bit_count() <= 3
               for v, row in enumerate(rows))
    assert all(bool(rows[v] >> w & 1) == bool(rows[w] >> v & 1)
               for v in range(n) for w in range(n))
    sources = independent_sources(rows)
    assert sources == [(1 << n) - 1], ('not strongly connected', instance.edges)
    colors, queue = {0: 0}, [0]
    for v in queue:
        for w in range(n):
            if rows[v] >> w & 1:
                if w not in colors:
                    colors[w] = 1 - colors[v]
                    queue.append(w)
                assert colors[v] != colors[w], ('not bipartite', instance.edges)
    assert n - match(rows, (1 << n) - 1) == 2 * len(instance.edges)


def girth(rows):
    """Undirected girth; return None for a forest."""
    best = None
    for start in range(len(rows)):
        distance, parent, queue = {start: 0}, {start: -1}, [start]
        for u in queue:
            for v in range(len(rows)):
                if not rows[u] >> v & 1:
                    continue
                if v not in distance:
                    distance[v] = distance[u] + 1
                    parent[v] = u
                    queue.append(v)
                elif v != parent[u]:
                    length = distance[u] + distance[v] + 1
                    best = length if best is None else min(best, length)
    return best


def check_source(edges, exhaustive_extras=False):
    instance = vertex_cover_graph(edges)
    structural_checks(instance)
    rows, n = instance.rows, len(instance.rows)
    forced = {v for pair in instance.forced_pairs for v in pair}
    assert sc(rows, forced)
    # Compute every one-state repair against the actual full state graph.
    for pair, expected in zip(instance.forced_pairs, instance.repair_sets):
        for a in pair:
            actual = {v for v in range(n) if v not in forced
                      and sc(rows, (forced - {a}) | {v})}
            assert actual == set(expected), (edges, pair, actual, expected)
    covers = source_covers(edges)
    optimum = min(map(len, covers))
    for cover in covers:
        selected = lift_cover(instance, cover)
        assert robust(rows, selected), (edges, cover, 'lift')
        decoded = extract_cover(instance, selected)
        assert len(decoded) <= len(selected) - 2 * len(edges)
        assert all(set(e) & decoded for e in edges)
    candidates = [v for v in range(n) if v not in forced and
                  v not in instance.centers]
    # Centers cannot repair any twin, as independently checked above.
    sizes = range(len(candidates) + 1) if exhaustive_extras else range(optimum + 1)
    target_optimum = None
    for size in sizes:
        for extra in combinations(candidates, size):
            selected = forced | set(extra)
            predicted = all(set(extra) & r for r in instance.repair_sets)
            actual = robust(rows, selected)
            assert actual == predicted, (edges, extra, 'full failure oracle')
            if actual:
                decoded = extract_cover(instance, selected)
                assert all(set(e) & decoded for e in edges)
                assert len(decoded) <= size
                target_optimum = min(target_optimum or size, size)
        if target_optimum is not None and not exhaustive_extras:
            break
    assert target_optimum == optimum, (edges, target_optimum, optimum)


def rank_table(instance):
    """All dual ranks on target states, using exact subset-DP matching ranks."""
    targets = [v for v in range(len(instance.rows)) if v not in instance.centers]
    index = {v: j for j, v in enumerate(targets)}
    rows = [sum(1 << index[v] for v in targets if instance.rows[u] >> v & 1)
            for u in sorted(instance.centers)]
    reachable = matching_subsets(rows)
    full = (1 << len(targets)) - 1
    ranks = [mask.bit_count() if mask in reachable else 0
             for mask in range(full + 1)]
    for j in range(len(targets)):
        for mask in range(full + 1):
            if mask >> j & 1:
                ranks[mask] = max(ranks[mask], ranks[mask ^ (1 << j)])
    nu = ranks[full]
    return targets, [s.bit_count() + ranks[full ^ s] - nu
                     for s in range(full + 1)]


def check_parallel_and_girth(edges, rng):
    base = vertex_cover_graph(edges)
    enlarged = subdivide_ports(base)
    structural_checks(enlarged)
    before, after = girth(base.rows), girth(enlarged.rows)
    assert after == (None if before is None else 3 * before)
    targets, ranks = rank_table(base)
    d = 2 * len(edges)
    distinct_opt = min(s.bit_count() for s, rank in enumerate(ranks) if rank == d
                       and all(ranks[s ^ (1 << j)] == d
                               for j in range(len(targets)) if s >> j & 1))
    # With repeated actuators, duplicate exactly the coloops of the support.
    repeated_opt = min(s.bit_count() + sum(ranks[s ^ (1 << j)] < d
                                         for j in range(len(targets)) if s >> j & 1)
                       for s, rank in enumerate(ranks) if rank == d)
    assert distinct_opt == repeated_opt == d + min(map(len, source_covers(edges)))
    full = (1 << len(enlarged.rows)) - 1
    nu = match(enlarged.rows, full)
    index = {v: j for j, v in enumerate(targets)}
    for _ in range(250):
        selected = {v for v in range(len(enlarged.rows)) if rng.randrange(2)}
        mask = sum(1 << v for v in selected)
        projected = {enlarged.projection[v] for v in selected}
        p = sum(1 << index[v] for v in projected if v in index)
        rank = len(selected) + match(enlarged.rows, full ^ mask) - nu
        assert rank == ranks[p], ('parallel-extension rank', edges, selected)
    best = min(source_covers(edges), key=len)
    assert robust(enlarged.rows, lift_cover(base, best))
    if len(edges) == 3 and len(source_covers(edges)) == 4:  # triangle
        twice = subdivide_ports(enlarged)
        assert girth(twice.rows) == 9 * girth(base.rows)
        assert robust(twice.rows, lift_cover(base, best))
        twice_full = (1 << len(twice.rows)) - 1
        twice_nu = match(twice.rows, twice_full)
        for _ in range(25):
            selected = {v for v in range(len(twice.rows)) if rng.randrange(2)}
            mask = sum(1 << v for v in selected)
            projected = {twice.projection[v] for v in selected}
            p = sum(1 << index[v] for v in projected if v in index)
            rank = len(selected) + match(twice.rows, twice_full ^ mask) - twice_nu
            assert rank == ranks[p]
    return len(ranks)


def compressed_triangle():
    # One twin per edge was deleted from the original, unported construction.
    edges = ((0, 1), (0, 2), (1, 2))
    rows = [0] * 15

    def link(u, v):
        rows[u] |= 1 << v
        rows[v] |= 1 << u

    for j, (u, v) in enumerate(edges):
        c, a, x = 3 * j, 3 * j + 1, 3 * j + 2
        link(c, a)
        link(c, x)
        for vertex in (u, v):
            link(9 + 2 * vertex, x)
    for v in range(3):
        link(9 + 2 * v, 10 + 2 * v)
    _, placements = brute_robust(tuple(rows))
    assert min(map(int.bit_count, placements)) == 4
    assert robust(rows, {1, 10, 12, 14})
    return len(rows), 1 << len(rows)


def check_leaf_completion():
    rng, count = Random(15513), 0
    for edges in (((0, 1),), ((0, 1), (1, 2)), ((0, 1), (2, 3)),
                  ((0, 1), (0, 2), (1, 2))):
        base = vertex_cover_graph(edges)
        completed = complete_leaves(base)
        padded = subdivide_ports(completed)
        assert girth(padded.rows) == 3 * girth(completed.rows)
        best = min(source_covers(edges), key=len)
        for graph in (completed, padded):
            structural_checks(graph)
            assert min(row.bit_count() for row in graph.rows) == 2
            assert robust(graph.rows, lift_cover(base, best))
            for _ in range(60):
                selected = {v for v in range(len(graph.rows)) if rng.randrange(2)}
                if not robust(graph.rows, selected):
                    continue
                normalized = normalize_subdivision(base, graph, selected)
                assert len(normalized) <= len(selected) and robust(base.rows, normalized)
                cover = extract_cover(base, normalized)
                assert len(cover) <= len(selected) - 2 * len(edges)
                count += 1
    print('Leafless cover gadgets: eight realizations,', count,
          'robust reverse maps PASS', flush=True)


if __name__ == '__main__':
    count = 0
    for n in range(2, 5):
        possible = list(combinations(range(n), 2))
        for bits in range(1, 1 << len(possible)):
            edges = tuple(e for j, e in enumerate(possible) if bits >> j & 1)
            check_source(edges, exhaustive_extras=(len(edges) <= 2))
            count += 1
        print('Source graphs through n=%d: %d PASS' % (n, count), flush=True)
    for edges in (tuple(combinations(range(5), 2)),
                  tuple((0, v) for v in range(1, 5)),
                  tuple((0, v) for v in range(1, 7))):
        check_source(edges)
    print('K5 and stars with 4 and 6 leaves: PASS', flush=True)
    for edges in (((0, 1),), ((0, 1), (1, 2))):
        instance = vertex_cover_graph(edges)
        _, placements = brute_robust(instance.rows)
        optimum = 2 * len(edges) + min(map(len, source_covers(edges)))
        assert min(map(int.bit_count, placements)) == optimum
        for s in placements:
            selected = {v for v in range(len(instance.rows)) if s >> v & 1}
            decoded = extract_cover(instance, selected)
            assert all(set(e) & decoded for e in edges)
    print('Unrestricted actuator enumeration: 2^10 + 2^18 placements PASS', flush=True)
    rng = Random(20261104)
    masks = 0
    for edges in (((0, 1),), ((0, 1), (1, 2)), ((0, 1), (0, 2), (1, 2)),
                  ((0, 1), (0, 2), (0, 3))):
        masks += check_parallel_and_girth(edges, rng)
    print('Parallel-extension audit: %d support subsets, 1025 full-state rank queries PASS'
          % masks, flush=True)
    vertices, subsets = compressed_triangle()
    print('Compressed triangle: n=%d, all %d placements, optimum 4 (not 5) PASS'
          % (vertices, subsets))
    check_leaf_completion()
