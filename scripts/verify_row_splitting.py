"""Exact rank, optimum and reverse-map checks for arbitrary row splitting."""
from itertools import combinations, combinations_with_replacement
from random import Random
from audit_hardness import augmented_matchable
from bounded_core_solver import matching_number, source_components
from hardness_gadgets import complete_leaves, subdivide_ports
from row_splitting import (cotransversal_from_placement, cotransversal_reduction,
                           avoid_presentation_element,
                           dual_rank_oracle, full_rank_presentation, normalize_row_splits,
                           ternary_presentation)
from verify_gammoid_obstructions import dual_ranks, robust, transversal_ranks


def audit_compression():
    count = 0
    for n in range(1, 4):
        for r in range(7):
            for rows in combinations_with_replacement(range(1 << n), r):
                original = transversal_ranks(rows, n)
                compressed = full_rank_presentation(rows, n)
                assert len(compressed) == original[-1]
                assert transversal_ranks(compressed, n) == original
                assert sum(row.bit_count() for row in compressed) <= sum(
                    row.bit_count() for row in rows)
                count += 1
    print('Full-rank presentation compression:', count, 'presentations PASS', flush=True)


def audit_local_deletion():
    sources = exchanges = 0
    for n in range(2, 5):
        full = (1 << n) - 1
        for r in range(1, n):
            for rows in combinations_with_replacement(range(1, 1 << n), r):
                tr = transversal_ranks(rows, n)
                if tr[-1] != r or any(tr[1 << v] == 0 for v in range(n)):
                    continue
                ranks = dual_ranks(rows, n)
                robust_masks = [mask for mask in range(1 << n) if robust(ranks, mask)]
                for x in range(n):
                    if sum(bool(row >> x & 1) for row in rows) < 2:
                        continue
                    assert min(mask.bit_count() for mask in robust_masks) == min(
                        mask.bit_count() for mask in robust_masks if not mask >> x & 1)
                    for mask in robust_masks:
                        selected = {v for v in range(n) if mask >> v & 1}
                        normalized = avoid_presentation_element(rows, n, selected, x)
                        assert x not in normalized and len(normalized) <= len(selected)
                        assert robust(ranks, sum(1 << v for v in normalized))
                        exchanges += 1
                sources += 1
    print('Individual column deletion:', sources, 'sources,', exchanges,
          'exchanges PASS', flush=True)


def audit():
    rng = Random(70462)
    count = reverse = restrictions = state_checks = 0

    def check(rows, n, exhaustive=False, graph_check=False, pad=False):
        nonlocal count, reverse, restrictions, state_checks
        rows = tuple(rows)
        ranks = dual_ranks(rows, n)
        split, size, steps = ternary_presentation(rows, n)
        d = n - len(rows)
        assert size - len(split) == d
        assert max(row.bit_count() for row in split) <= 3
        assert len(steps) == sum(max(0, row.bit_count() - 3) for row in rows)
        final = dual_ranks(split, size) if exhaustive else None
        oracle = dual_rank_oracle(split, size)
        for mask, r in enumerate(ranks):
            selected = {v for v in range(n) if mask >> v & 1}
            assert oracle(selected) == r
            restrictions += 1
        tr = transversal_ranks(rows, n)
        if all(tr[(1 << u) | (1 << v)] == 2 for u, v in combinations(range(n), 2)):
            new_tr = transversal_ranks(split, size)
            assert all(new_tr[(1 << u) | (1 << v)] == 2
                       for u, v in combinations(range(size), 2))
        supports = list(range(1 << size)) if exhaustive else []
        supports += [sum(1 << v for v in range(size) if rng.randrange(2))
                     for _ in range(40)]
        for mask in supports:
            selected = {v for v in range(size) if mask >> v & 1}
            valid = (robust(final, mask) if final is not None else
                     oracle(selected) == d and all(oracle(selected - {v}) == d
                                                   for v in selected))
            if not valid:
                continue
            normalized = normalize_row_splits(rows, n, selected)
            assert len(normalized) <= len(selected)
            assert robust(ranks, sum(1 << v for v in normalized))
            reverse += 1
        if exhaustive:
            assert min(mask.bit_count() for mask in range(1 << n) if robust(ranks, mask)) == min(
                mask.bit_count() for mask in range(1 << size) if robust(final, mask))
        if graph_check:
            graph = cotransversal_reduction(rows, n)
            if pad:
                graph = subdivide_ports(graph)
            total = len(graph.rows)
            full = (1 << total) - 1
            nu = matching_number(graph.rows, full)
            assert total - nu == d
            assert max(row.bit_count() for row in graph.rows) <= 3
            assert len(source_components(graph.rows)) == 1
            for mask in range(1 << n):
                if not robust(ranks, mask):
                    continue
                selected = {v for v in range(n) if mask >> v & 1}
                assert all(augmented_matchable(graph.rows, selected - {v}) for v in selected)
                state_checks += len(selected)
            for _ in range(10):
                selected = {v for v in range(total) if rng.randrange(2)}
                if not selected or not all(augmented_matchable(graph.rows, selected - {v})
                                           for v in selected):
                    continue
                normalized = cotransversal_from_placement(rows, n, graph, selected)
                assert len(normalized) <= len(selected)
                assert robust(ranks, sum(1 << v for v in normalized))
                reverse += 1
        count += 1

    for n in range(2, 5):
        full = (1 << n) - 1
        for r in range(1, n):
            for rows in combinations_with_replacement(range(1, 1 << n), r):
                if max(row.bit_count() for row in rows) <= 3:
                    continue
                covered = 0
                for row in rows:
                    covered |= row
                if covered == full and transversal_ranks(rows, n)[-1] == r:
                    check(rows, n, exhaustive=True)
        print('Arbitrary row splitting through', n, 'elements:', count, 'sources PASS', flush=True)
    # Contract degree-two columns of the certified nontransversal rank-five
    # prototype: merging their two rows produces genuine size-four sources.
    prototype = (1033, 4612, 2084, 1064, 560, 352, 448, 4114)
    for x in (2, 3, 4, 6, 8, 9, 10, 12):
        incident = [j for j, row in enumerate(prototype) if row >> x & 1]
        combined = 0
        rows = []
        for j, row in enumerate(prototype):
            if j in incident:
                combined |= row
            else:
                rows.append(row)
        rows.append(combined ^ (1 << x))
        rows = tuple((row & ((1 << x) - 1)) | ((row >> (x + 1)) << x) for row in rows)
        ranks = dual_ranks(rows, 12)
        assert min(mask.bit_count() for mask in range(1 << 12) if robust(ranks, mask)) == 7
        check(rows, 12)
    for n in range(5, 9):
        for _ in range(40):
            r = rng.randrange(1, min(n, 5))
            # One full row forces a nontrivial split and covers every element.
            rows = ((1 << n) - 1,) + tuple(rng.randrange(1, 1 << n) for _ in range(r - 1))
            if transversal_ranks(rows, n)[-1] != r:
                continue
            split, size, _ = ternary_presentation(rows, n)
            if size > 17:
                continue
            check(rows, n, exhaustive=size <= 10)
    for rows, n, pad in (((15,), 4, False), ((15, 240), 8, True),
                         ((31, 3, 12), 5, False)):
        check(rows, n, exhaustive=True, graph_check=True, pad=pad)
    print('Arbitrary row splitting:', count, 'sources,', restrictions,
          'restriction ranks,', reverse, 'reverse maps,', state_checks,
          'direct input-failure checks PASS', flush=True)


def audit_control_sources():
    # These are actual directed source graphs with redundant presentation rows.
    sources = [(14, 1, 1, 1), (30, 1, 1, 1, 1), (30, 3, 1, 1, 1),
               (6, 8, 8, 1), (14, 5, 1, 1)]
    count = 0
    rng = Random(513)
    for rows in sources:
        n = len(rows)
        assert source_components(rows) == [(1 << n) - 1]
        nu = matching_number(rows, (1 << n) - 1)
        if nu == n:
            continue  # The zero-deficiency case has optimum two directly.
        rank = dual_rank_oracle(full_rank_presentation(rows, n), n)
        source_opt = min(mask.bit_count() for mask in range(1 << n)
                         if len(selected := {v for v in range(n) if mask >> v & 1}) >= 2
                         and rank(selected) == n - nu
                         and all(rank(selected - {v}) == n - nu for v in selected))
        graph = cotransversal_reduction(rows, n)
        for mask in range(1 << n):
            selected = {v for v in range(n) if mask >> v & 1}
            if len(selected) < 2 or not all(augmented_matchable(rows, selected - {v})
                                            for v in selected):
                continue
            assert all(augmented_matchable(graph.rows, selected - {v}) for v in selected)
        for _ in range(20):
            selected = {v for v in range(len(graph.rows)) if rng.randrange(2)}
            if not selected or not all(augmented_matchable(graph.rows, selected - {v})
                                       for v in selected):
                continue
            normalized = cotransversal_from_placement(rows, n, graph, selected)
            assert source_opt <= len(normalized) <= len(selected)
            assert all(augmented_matchable(rows, normalized - {v}) for v in normalized)
        count += 1
    print('Actual strongly connected directed sources:', count, 'PASS', flush=True)


def audit_leaf_completion():
    rng = Random(1265)
    count = queries = reverse = 0
    fixtures = (((3,), 2), ((7,), 3), ((3, 12), 4), ((3, 4), 3),
                ((15,), 4), ((31, 3, 12), 5),
                ((1033, 4612, 2084, 1064, 560, 352, 448, 4114), 13))
    for rows, n in fixtures:
        for pad in (False, True):
            graph = complete_leaves(cotransversal_reduction(rows, n))
            if pad:
                graph = subdivide_ports(graph)
            total, full = len(graph.rows), (1 << len(graph.rows)) - 1
            assert all(2 <= row.bit_count() <= 3 for row in graph.rows)
            assert source_components(graph.rows) == [full]
            nu = matching_number(graph.rows, full)
            rank = dual_rank_oracle(rows, n)
            split, size, _ = ternary_presentation(rows, n)
            split_rank = dual_rank_oracle(split, size)
            assert total - nu == n - len(rows)
            for _ in range(80):
                selected = {v for v in range(total) if rng.randrange(2)}
                support = {graph.projection[v] for v in selected if graph.projection[v] >= 0}
                mask = sum(1 << v for v in selected)
                actual = len(selected) + matching_number(graph.rows, full ^ mask) - nu
                assert actual == split_rank(support)
                queries += 1
                if (actual != n - len(rows) or any(
                        len(selected) - 1 + matching_number(graph.rows, full ^ (mask ^ (1 << v)))
                        - nu != actual for v in selected)):
                    continue
                original = cotransversal_from_placement(rows, n, graph, selected)
                assert len(original) <= len(selected)
                assert rank(original) == actual and all(rank(original - {v}) == actual for v in original)
                reverse += 1
            # Explicit failures after deliberately omitting a prescribed state
            # certify absence of mandatory actuators in the actual control model.
            for excluded in range(min(5, total)):
                selected = set(range(total)) - {excluded}
                for failed in sorted(selected)[:3]:
                    assert augmented_matchable(graph.rows, selected - {failed})
            count += 1
    print('Leaf completion:', count, 'realizations,', queries,
          'full-state ranks,', reverse, 'reverse maps PASS', flush=True)


if __name__ == '__main__':
    audit_compression()
    audit_local_deletion()
    audit()
    audit_control_sources()
    audit_leaf_completion()
