"""Exact finite audits for deterministic actuator placement and complexity barriers."""

from collections import deque
from fractions import Fraction
from itertools import combinations, product
from math import ceil, comb
from random import Random

from audit_linear_accessibility_quotas import graph_survives, optimum, small_rank, subsets
from audit_weighted_quotas import core as mixed_core
from audit_mixed_accessibility_quotas import group_prices
from bounded_core_solver import matching_number, source_components
from deterministic_actuator_placement import (
    mimicking_digraph,
    rank_two_representation,
    solve_deficiency_one,
    solve_deficiency_two,
    solve_deficiency_zero,
    solve_small_deficiency,
    sparse_paving_lower_bound_instance,
    strict_gammoid_digraph,
    strict_gammoid_q_parity_to_actuator,
)
from randomized_actuator_placement import verify_placement


def rect_matching_number(rows, n_right, allowed):
    """Exact matching size for a rectangular bipartite graph with n_right columns."""
    pad = n_right - len(rows)
    return matching_number(rows + (0,) * pad if pad > 0 else rows, allowed)


def dual_matching_rank(rows, subset):
    """Exact rank of subset in M = T(B_G)^* via bipartite matching."""
    n = len(rows)
    full = (1 << n) - 1
    rho = matching_number(rows, full)
    mask = sum(1 << e for e in subset)
    return len(subset) - rho + matching_number(rows, full ^ mask)


def max_vertex_disjoint_paths(adj, sources, sinks):
    """Independent max-flow vertex-disjoint path count from sources to sinks."""
    sources = frozenset(sources)
    sinks = frozenset(sinks)
    if not sources or not sinks:
        return 0
    n = len(adj)
    src = 2 * n
    snk = 2 * n + 1
    total = 2 * n + 2
    cap = [dict() for _ in range(total)]

    def add_edge(u, v, c):
        cap[u][v] = cap[u].get(v, 0) + c
        cap[v].setdefault(u, 0)

    for v in range(n):
        add_edge(2 * v, 2 * v + 1, 1)
    for u, mask in enumerate(adj):
        while mask:
            bit = mask & -mask
            mask ^= bit
            v = bit.bit_length() - 1
            add_edge(2 * u + 1, 2 * v, 1)
    for s in sources:
        add_edge(src, 2 * s, 1)
    for t in sinks:
        add_edge(2 * t + 1, snk, 1)

    flow = 0
    while True:
        parent = [-1] * total
        parent[src] = src
        q = deque([src])
        while q and parent[snk] == -1:
            u = q.popleft()
            for v, c in cap[u].items():
                if c > 0 and parent[v] == -1:
                    parent[v] = u
                    q.append(v)
        if parent[snk] == -1:
            break
        v = snk
        while v != src:
            u = parent[v]
            cap[u][v] -= 1
            cap[v][u] += 1
            v = u
        flow += 1
    return flow


def min_ab_vertex_cut(adj_sets, a_set, b_set, preferred_set):
    """Minimum-cardinality (A, B)-vertex cut maximizing overlap with preferred_set."""
    n = len(adj_sets)
    if not a_set or not b_set:
        return frozenset()
    src = 2 * n
    snk = 2 * n + 1
    total = 2 * n + 2
    big = n + 1
    cap = [dict() for _ in range(total)]

    def add_edge(u, v, c):
        cap[u][v] = cap[u].get(v, 0) + c
        cap[v].setdefault(u, 0)

    for v in range(n):
        weight = n if v in preferred_set else big
        add_edge(2 * v, 2 * v + 1, weight)
    inf_cap = (n + 2) * big
    for u in range(n):
        for v in adj_sets[u]:
            add_edge(2 * u + 1, 2 * v, inf_cap)
    for a in a_set:
        add_edge(src, 2 * a, inf_cap)
    for b in b_set:
        add_edge(2 * b + 1, snk, inf_cap)

    while True:
        parent = [-1] * total
        parent[src] = src
        q = deque([src])
        while q and parent[snk] == -1:
            u = q.popleft()
            for v, c in cap[u].items():
                if c > 0 and parent[v] == -1:
                    parent[v] = u
                    q.append(v)
        if parent[snk] == -1:
            break
        v = snk
        bottleneck = inf_cap
        while v != src:
            u = parent[v]
            bottleneck = min(bottleneck, cap[u][v])
            v = u
        v = snk
        while v != src:
            u = parent[v]
            cap[u][v] -= bottleneck
            cap[v][u] += bottleneck
            v = u

    reachable = {src}
    q = deque([src])
    while q:
        u = q.popleft()
        for v, c in cap[u].items():
            if c > 0 and v not in reachable:
                reachable.add(v)
                q.append(v)

    cut = {v for v in range(n) if (2 * v) in reachable and (2 * v + 1) not in reachable}
    return frozenset(cut)


def exact_cut_covering_set(adj, terminal_subset, sinks):
    """Compute a small exact cut-covering set Z for (S, T) in D~_H."""
    n = len(adj)
    terminal_list = sorted(terminal_subset)
    d = len(sinks)
    s_indices = tuple(range(n, n + len(terminal_list)))
    t_indices = tuple(range(n + len(terminal_list), n + len(terminal_list) + d))
    total = n + len(terminal_list) + d

    adj_sets = [set() for _ in range(total)]
    for u, mask in enumerate(adj):
        while mask:
            bit = mask & -mask
            mask ^= bit
            adj_sets[u].add(bit.bit_length() - 1)
    for pos, w in enumerate(terminal_list):
        adj_sets[n + pos].add(w)
    for pos, t in enumerate(sinks):
        adj_sets[t].add(n + len(terminal_list) + pos)

    z = set(s_indices) | set(t_indices)
    pairs = [
        (a, b)
        for sz_a in range(1, len(s_indices) + 1)
        for a in combinations(s_indices, sz_a)
        for sz_b in range(1, len(t_indices) + 1)
        for b in combinations(t_indices, sz_b)
    ]
    # Process larger terminal pairs first to reuse cuts.
    pairs.sort(key=lambda ab: -(len(ab[0]) + len(ab[1])))
    for a, b in pairs:
        cut = min_ab_vertex_cut(adj_sets, a, b, z)
        z.update(cut)
    return frozenset(z)


def audit_small_deficiency():
    rng = Random(20261126)
    checked = 0
    counts_by_d = {0: 0, 1: 0, 2: 0}

    # Explicit adversarial cases covering d in {0, 1, 2}.
    adversarial = [
        # d=0, single state without self-loop, k=0 and k=1.
        ((0,), (5,), 0, {0}),
        ((0,), (5,), 1, {0}),
        ((0,), (5,), 1, set()),
        # d=1, multiple source SCCs, loops and nonloops, safe and vulnerable.
        ((2, 1, 0, 12, 8), (3, 4, 2, 1, 5), 1, {0, 2, 3}),
        ((2, 1, 0, 12, 8), (3, 4, 2, 1, 5), 2, {0, 1, 2, 3, 4}),
        # d=2, parallel classes and loops with zero weights.
        ((0b01000, 0b01000, 0b10000, 0b10000, 0b00001), (0, 2, 1, 4, 4), 1, {0, 1, 2}),
        ((0b00100, 0b01000, 0b10000, 0b00100, 0b01000), (1, 2, 3, 0, 0), 2, {0, 1, 2, 3, 4}),
    ]
    for rows, weights, k, vuln in adversarial:
        n = len(rows)
        d = n - matching_number(rows, (1 << n) - 1)
        expected = optimum(n, weights, k, vuln, lambda ch: graph_survives(rows, ch))
        sol = solve_small_deficiency(rows, weights, k, vuln)
        cost = None if sol is None else sum(weights[e] for e in sol.placement)
        assert cost == expected
        if sol is not None:
            assert verify_placement(rows, k, vuln, sol.placement, sol.witness)
        if d == 2:
            cols = rank_two_representation(rows)
            for sub in subsets(range(n)):
                assert small_rank(cols[e] for e in sub) == dual_matching_rank(rows, sub)
        counts_by_d[d] += 1
        checked += 1

    # Random exhaustive-comparison sweep across d in {0, 1, 2}.
    while min(counts_by_d.values()) < 85:
        n = rng.randrange(1, 7)
        rows = tuple(rng.randrange(1 << n) for _ in range(n))
        rho = matching_number(rows, (1 << n) - 1)
        d = n - rho
        if d not in (0, 1, 2):
            continue
        if counts_by_d[d] >= 95:
            continue
        weights = tuple(rng.randrange(7) for _ in range(n))
        vuln = frozenset(e for e in range(n) if rng.randrange(2))
        k = rng.randrange(3 if d == 2 else 4)
        expected = optimum(n, weights, k, vuln, lambda ch: graph_survives(rows, ch))
        sol = solve_small_deficiency(rows, weights, k, vuln)
        cost = None if sol is None else sum(weights[e] for e in sol.placement)
        assert cost == expected, (rows, weights, k, vuln, cost, expected)
        if sol is not None:
            assert verify_placement(rows, k, vuln, sol.placement, sol.witness)
            for fail in subsets(sol.placement & vuln, k):
                assert graph_survives(rows, sol.placement - fail)
        if d == 2:
            cols = rank_two_representation(rows)
            for sub in subsets(range(n)):
                assert small_rank(cols[e] for e in sub) == dual_matching_rank(rows, sub)
        counts_by_d[d] += 1
        checked += 1

    # Larger d = 1 graphs with multiple protected source SCCs (n = 18, u = 5, k = 3).
    large_d1_checked = 0
    for trial in range(15):
        n = 18
        rows = [0] * n
        # Five 3-cycles (states 0..14) as source SCCs plus three sink/intermediate states 15..17.
        for g in range(5):
            b = 3 * g
            rows[b] |= 1 << (b + 1)
            rows[b + 1] |= 1 << (b + 2)
            rows[b + 2] |= 1 << b
            # Arcs into intermediate chain 15 -> 16 -> 17 (leaving 17 unmatched).
            rows[b] |= 1 << 15
            if rng.randrange(2):
                rows[b + 1] |= 1 << 16
        rows[15] |= 1 << 16
        rows[16] |= 1 << 17
        rows = tuple(rows)
        assert n - matching_number(rows, (1 << n) - 1) == 1
        weights = tuple(rng.randrange(1, 20) for _ in range(n))
        vuln = frozenset(e for e in range(n) if e % 3 != 0 or rng.randrange(2))
        k = 3
        sol = solve_deficiency_one(rows, weights, k, vuln)
        if sol is not None:
            assert verify_placement(rows, k, vuln, sol.placement, sol.witness)
            # Compare against explicit group-choice enumeration.
            adj, sinks, _ = strict_gammoid_digraph(rows)
            nonloops = {
                e for e in range(n)
                if max_vertex_disjoint_paths(adj, (e,), sinks) == 1
            }
            sources = [
                tuple(e for e in range(n) if mask >> e & 1)
                for mask in source_components(rows)
            ]
            # Check each source SCC admits a valid option and compare cost on feasible small subsets.
            opt_check = _d1_group_bruteforce(sources, set(range(n)) - set().union(*sources),
                                             weights, vuln, nonloops, k + 1)
            assert sum(weights[e] for e in sol.placement) == opt_check
        large_d1_checked += 1

    return checked, counts_by_d, large_d1_checked


def _d1_group_bruteforce(groups, ungrouped, weights, vuln, nonloops, threshold):
    group_opts = []
    for g in groups:
        opts = []
        for sub in subsets(g):
            safe_sub = sub - vuln
            vuln_sub = sub & vuln
            if safe_sub or len(vuln_sub) >= threshold:
                opts.append(sub)
        if not opts:
            return None
        group_opts.append(opts)
    ung_opts = list(subsets(ungrouped))
    best = None
    for prod in product(*group_opts, ung_opts):
        chosen = frozenset().union(*prod)
        if (chosen - vuln) & nonloops or len((chosen & vuln) & nonloops) >= threshold:
            cost = sum(weights[e] for e in chosen)
            if best is None or cost < best:
                best = cost
    return best


def _f_p_rank(matrix, p):
    rows = [row[:] for row in matrix]
    if not rows or not rows[0]:
        return 0
    m = len(rows)
    ncols = len(rows[0])
    r = 0
    for col in range(ncols):
        piv = next((i for i in range(r, m) if rows[i][col] % p != 0), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][col], -1, p)
        prow = [(val * inv) % p for val in rows[r][col:]]
        for i in range(r + 1, m):
            factor = rows[i][col] % p
            if factor:
                rows[i][col:] = [(a - factor * b) % p for a, b in zip(rows[i][col:], prow)]
        r += 1
        if r == m:
            break
    return r


def _f_p_nullspace(matrix, p):
    rows = [[v % p for v in row] for row in matrix]
    m = len(rows)
    ncols = len(rows[0]) if m else 0
    pivots = []
    r = 0
    for col in range(ncols):
        piv = next((i for i in range(r, m) if rows[i][col] != 0), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][col], -1, p)
        rows[r][col:] = [(val * inv) % p for val in rows[r][col:]]
        prow = rows[r][col:]
        for i in range(m):
            if i != r and rows[i][col] != 0:
                factor = rows[i][col]
                rows[i][col:] = [(a - factor * b) % p for a, b in zip(rows[i][col:], prow)]
        pivots.append(col)
        r += 1
    free_cols = [c for c in range(ncols) if c not in pivots]
    basis_rows = []
    for fc in free_cols:
        vec = [0] * ncols
        vec[fc] = 1
        for i, pc in enumerate(pivots):
            vec[pc] = (-rows[i][fc]) % p
        basis_rows.append(vec)
    return basis_rows


def _is_prime(x):
    if x < 2:
        return False
    d = 2
    while d * d <= x:
        if x % d == 0:
            return False
        d += 1
    return True


def _next_prime(x):
    cand = max(2, x + 1)
    while not _is_prime(cand):
        cand += 1
    return cand


def _koana_kumabe_check(b_rows, n_right, delta_num, delta_den, *, use_rho_scale=False):
    """Exact F_p check of Koana-Kumabe transversal and dual almost-multilinear matrix."""
    n_left = len(b_rows)
    n_total = n_left + n_right
    rho = rect_matching_number(b_rows, n_right, (1 << n_right) - 1)
    scale_factor = rho if use_rho_scale else n_total
    t = ceil(2 * n_total * scale_factor * delta_den / delta_num)
    ell = t + 2 * n_total
    big_k = n_total * (t + n_total) + 1
    p = _next_prime(big_k + t - 1)
    alphas = list(range(1, n_right + 1))
    betas = list(range(n_right + 1, n_total + 1))

    # Build matrix C of size (ell * n_left) x (t * n_right) over F_p.
    # For each edge (w, v), expand (beta_w - alpha_v + Z)^K * (beta_w + Z)^a mod Z^ell.
    nrows = ell * n_left
    ncols = t * n_right
    c_mat = [[0] * ncols for _ in range(nrows)]

    for w in range(n_left):
        beta = betas[w]
        for v in range(n_right):
            if not (b_rows[w] >> v & 1):
                continue
            alpha = alphas[v]
            # Compute (Z + (beta - alpha))^K mod Z^ell via binomial coefficients.
            diff = (beta - alpha) % p
            # Start with poly_k[h] = binom(K, h) * diff^(K-h) mod p.
            comb_k = 1
            poly_k = [0] * ell
            for h in range(min(ell, big_k + 1)):
                poly_k[h] = (comb_k * pow(diff, big_k - h, p)) % p
                comb_k = (comb_k * (big_k - h) * pow(h + 1, -1, p)) % p
            # Multiply sequentially by (beta + Z) for a = 0 .. t - 1.
            cur = poly_k
            for a in range(t):
                row_base = w * ell
                col_idx = v * t + a
                for h in range(ell):
                    c_mat[row_base + h][col_idx] = cur[h]
                if a + 1 < t:
                    nxt = [0] * ell
                    for h in range(ell):
                        val = (beta * cur[h]) % p
                        if h > 0:
                            val = (val + cur[h - 1]) % p
                        nxt[h] = val
                    cur = nxt

    # Check Claim 4.1 and Claim 4.2 on all subsets F of V_right.
    delta = Fraction(delta_num, delta_den)
    for sub in subsets(range(n_right)):
        r_t = rect_matching_number(b_rows, n_right, sum(1 << e for e in sub))
        cols_idx = [v * t + a for v in sub for a in range(t)]
        sub_mat = [[row[c] for c in cols_idx] for row in c_mat] if cols_idx else []
        rk = _f_p_rank(sub_mat, p) if cols_idx else 0
        assert t * r_t <= rk <= (t + 2 * n_total) * r_t
        assert abs(Fraction(rk, t) - r_t) <= delta

    # Dualize via Lemma 4.5 and verify 2*delta almost-multilinearity for T(B)^*.
    h_mat = _f_p_nullspace(c_mat, p)
    rho = rect_matching_number(b_rows, n_right, (1 << n_right) - 1)
    d_rank = n_right - rho
    for sub in subsets(range(n_right)):
        r_m = len(sub) - rho + rect_matching_number(
            b_rows, n_right, ((1 << n_right) - 1) ^ sum(1 << e for e in sub)
        )
        cols_idx = [v * t + a for v in sub for a in range(t)]
        sub_h = [[row[c] for c in cols_idx] for row in h_mat] if cols_idx else []
        rk_h = _f_p_rank(sub_h, p) if cols_idx else 0
        assert abs(Fraction(rk_h, t) - r_m) <= 2 * delta
        if r_m == d_rank:
            assert Fraction(rk_h, t) >= d_rank - 2 * delta
        else:
            assert Fraction(rk_h, t) <= (d_rank - 1) + 2 * delta
    return p, t, len(c_mat), len(c_mat[0])


def audit_gammoid_and_cut_covering():
    rng = Random(20261127)
    ip_checked = 0
    for _ in range(65):
        n = rng.randrange(4, 9)
        rows = tuple(rng.randrange(1 << n) for _ in range(n))
        adj, sinks, _ = strict_gammoid_digraph(rows)
        for sub in subsets(range(n)):
            expected = dual_matching_rank(rows, sub)
            got = max_vertex_disjoint_paths(adj, sub, sinks)
            assert got == expected
        ip_checked += 1

    # Mimicking digraph and flat-saturation core equivalence checks.
    mimic_checked = 0
    strict_reduction_count = 0
    for trial in range(30):
        # Construct state graphs with multi-hop routing paths and d in {1, 2, 3}.
        n = rng.randrange(7, 12)
        d = rng.randrange(1, min(4, n - 3))
        rho = n - d
        rows = [0] * n
        # First rho vertices are matched 1-to-1 to create an underlying DAG + cycles.
        for i in range(rho):
            rows[i] |= 1 << i
            # Forward arcs to higher indices + occasional back arcs.
            for j in range(n):
                if j != i and rng.random() < 0.3:
                    rows[i] |= 1 << j
        # Last d vertices are unmatched on the right if no row reaches them?
        # Wait: ensure left rows >= rho have empty rows so matching size <= rho is easy,
        # or add random arcs only from first rho rows so exact rho is guaranteed!
        rows = [0] * n
        for i in range(rho):
            rows[i] |= 1 << i
            for j in range(n):
                if j != i and rng.random() < 0.35:
                    rows[i] |= 1 << j
        rows = tuple(rows)
        assert matching_number(rows, (1 << n) - 1) == rho

        adj, sinks, _ = strict_gammoid_digraph(rows)
        # Compute a flat-saturation core W from mixed-accessibility-quotas on M = T^*.
        weights = [rng.randrange(1, 15) for _ in range(n)]
        colors = [rng.randrange(4) for _ in range(n)]
        renum = {g: idx for idx, g in enumerate(sorted(set(colors) | {0}))}
        colors = [renum[g] for g in colors]
        safe = [rng.randrange(3) == 0 for _ in range(n)]
        k = 1
        prices, protected = group_prices(weights, colors, safe, k + 1)
        if prices is None:
            # Make all elements safe in infeasible colors to ensure core construction runs.
            safe = [True] * n
            prices, protected = group_prices(weights, colors, safe, k + 1)

        w_core, _ = mixed_core(
            lambda sub: dual_matching_rank(rows, sub),
            d,
            k,
            weights,
            colors,
            safe,
            {},
            prices=prices,
            protected=protected,
        )
        if not w_core:
            continue
        # Test exact cut-covering set Z on subsets W of size <= min(4, len(w_core)).
        if len(w_core) <= 5:
            w_sub = frozenset(w_core)
        else:
            w_sub = frozenset(sorted(w_core)[:4])

        z_cover = exact_cut_covering_set(adj, w_sub, sinks)
        torso_adj, b_z_rows, idx_map = mimicking_digraph(adj, sinks, w_sub, z_cover)
        rho_z = len(b_z_rows)
        n_z = len(idx_map["z_ordered"])
        assert rho_z == n_z - d
        assert rect_matching_number(b_z_rows, n_z, (1 << n_z) - 1) == rho_z
        if rho_z < rho:
            strict_reduction_count += 1

        # Verify exact rank identity: T(B_Z)^*|_S = M|_W on all 2^|W| subsets!
        full_bz = (1 << n_z) - 1
        for sub in subsets(sorted(w_sub)):
            s_mask = sum(1 << idx_map["w_to_b_z"][w] for w in sub)
            r_bz = len(sub) - rho_z + rect_matching_number(b_z_rows, n_z, full_bz ^ s_mask)
            r_orig = dual_matching_rank(rows, sub)
            assert r_bz == r_orig

        # Verify optimum preservation when V \ W becomes loops in M'.
        if len(w_core) <= 5 and n <= 9:
            vuln = frozenset(e for e in range(n) if not safe[e])
            groups = [
                frozenset(e for e in range(n) if colors[e] == g)
                for g in set(colors) - {0}
            ]
            def orig_surv(chosen):
                return all(chosen & g for g in groups) and dual_matching_rank(rows, chosen) == d
            def compressed_surv(chosen):
                if not all(chosen & g for g in groups):
                    return False
                sub = chosen & w_sub
                s_mask = sum(1 << idx_map["w_to_b_z"][w] for w in sub)
                r_bz = len(sub) - rho_z + rect_matching_number(b_z_rows, n_z, full_bz ^ s_mask)
                return r_bz == d
            opt_orig = optimum(n, weights, k, vuln, orig_surv)
            opt_comp = optimum(n, weights, k, vuln, compressed_surv)
            assert opt_orig == opt_comp
        mimic_checked += 1

    assert strict_reduction_count > 0

    # Koana-Kumabe finite-field Wronskian construction & Lemma 4.5 verification.
    kk_checks = []
    # Case 1: |V|=3, |W|=1, rho=1, d=2, N_0=4, delta=1/5 (literal t = ceil(2*N_0^2/delta)).
    kk_checks.append(_koana_kumabe_check((0b110,), 3, 1, 5))
    # Case 2: |V|=3, |W|=2, rho=2, d=1, N_0=5, delta=1/5 (t = ceil(2*N_0*rho/delta)).
    kk_checks.append(_koana_kumabe_check((0b110, 0b011), 3, 1, 5, use_rho_scale=True))
    return ip_checked, mimic_checked, strict_reduction_count, kk_checks


def audit_sparse_paving_lower_bound():
    configs = [
        (3, 2, 1, False),
        (2, 2, 2, False),
        (3, 3, 1, False),
        (3, 4, 1, False),
        (3, 3, 1, True),
        (2, 3, 2, True),
        (3, 5, 1, True),
        (3, 6, 1, True),
    ]
    checked_instances = 0
    axiom_checks = 0
    indistinguishable_pairs = 0

    for u, r, k, loopless in configs:
        t = ((r - 1) // 2) if loopless else (r // 2)
        n, r_all, w, vuln, groups, opt_yes, opt_no, _ = sparse_paving_lower_bound_instance(
            u, r, k, None, loopless=loopless
        )
        assert n <= 12
        # Exhaustive local matroid rank axiom check across all 2^n subsets.
        rank_table = [r_all(frozenset(e for e in range(n) if m >> e & 1)) for m in range(1 << n)]
        assert rank_table[0] == 0
        assert rank_table[(1 << n) - 1] == r
        if loopless:
            assert all(rank_table[1 << e] == 1 for e in range(n))
        for m in range(1 << n):
            rm = rank_table[m]
            absent = [e for e in range(n) if not (m >> e & 1)]
            for e in absent:
                re = rank_table[m | (1 << e)]
                assert rm <= re <= rm + 1
                axiom_checks += 1
            for e, f in combinations(absent, 2):
                assert rank_table[m | (1 << e) | (1 << f)] + rm <= (
                    rank_table[m | (1 << e)] + rank_table[m | (1 << f)]
                )
                axiom_checks += 1

        def survives_one(survivor, tbl):
            if not all(survivor & g for g in groups):
                return False
            smask = sum(1 << e for e in survivor)
            return tbl[smask] == r

        got_no = optimum(n, w, k, vuln, lambda ch: survives_one(ch, rank_table))
        assert got_no == opt_no
        checked_instances += 1

        for hidden in combinations(range(u), t):
            n2, r_hid, w2, vuln2, groups2, opt_yes2, opt_no2, h_set = (
                sparse_paving_lower_bound_instance(u, r, k, hidden, loopless=loopless)
            )
            assert (n2, w2, vuln2, groups2, opt_yes2, opt_no2) == (
                n, w, vuln, groups, opt_yes, opt_no
            )
            hid_table = [
                r_hid(frozenset(e for e in range(n) if m >> e & 1)) for m in range(1 << n)
            ]
            assert hid_table[0] == 0 and hid_table[(1 << n) - 1] == r
            for m in range(1 << n):
                rm = hid_table[m]
                absent = [e for e in range(n) if not (m >> e & 1)]
                for e in absent:
                    re = hid_table[m | (1 << e)]
                    assert rm <= re <= rm + 1
                for e, f in combinations(absent, 2):
                    assert hid_table[m | (1 << e) | (1 << f)] + rm <= (
                        hid_table[m | (1 << e)] + hid_table[m | (1 << f)]
                    )
            # Check exact oracle indistinguishability away from H_{I*}.
            diff_count = sum(1 for m in range(1 << n) if hid_table[m] != rank_table[m])
            non_core_count = n - (2 * u + k)
            assert diff_count == (1 << non_core_count)
            got_yes = optimum(n, w, k, vuln, lambda ch: survives_one(ch, hid_table))
            assert got_yes == opt_yes
            indistinguishable_pairs += 1
            checked_instances += 1

    return checked_instances, axiom_checks, indistinguishable_pairs


def _solve_strict_gammoid_q_parity_exact(adj, sinks, groups):
    """Brute-force YES/NO oracle for Strict Gammoid q-Parity."""
    q = len(groups[0])
    d = len(sinks)
    t = d // q
    for chosen_groups in combinations(groups, t):
        union = tuple(v for g in chosen_groups for v in g)
        if max_vertex_disjoint_paths(adj, union, sinks) == d:
            return chosen_groups
    return None


def audit_strict_gammoid_q_parity_reduction():
    test_cases = [
        # 1. YES: d=3, t=1, g=2. G0={0,1,2} reaches {3,4,5}, G1={0,1,2} no -> disjoint!
        (6, [(0, 3), (1, 4), (2, 5)], (3, 4, 5), ((0, 1, 2), (3, 4, 5)), True),
        # 2. YES: d=3, t=1, g=2, G0=(0,1,2) blocked at 6, G1=(3,4,5) reaches {6,7,8}.
        (9, [(0, 6), (1, 6), (2, 7), (3, 6), (4, 7), (5, 8)], (6, 7, 8),
         ((0, 1, 2), (3, 4, 5)), True),
        # 3. NO: d=3, t=1, g=2, both groups route through a 2-vertex bottleneck {6, 7}.
        (11, [(i, 6) for i in range(6)] + [(j, 7) for j in range(6)]
         + [(6, 8), (6, 9), (7, 9), (7, 10)], (8, 9, 10),
         ((0, 1, 2), (3, 4, 5)), False),
        # 4. YES: d=3, t=1, g=3, cyclic digraph with 2-hop routing.
        (9, [(0, 3), (3, 0), (1, 4), (4, 6), (2, 5), (5, 7), (6, 7), (7, 8), (3, 6), (0, 8)],
         (6, 7, 8), ((0, 1, 2), (3, 4, 5)), True),
        # 5. NO: d=3, t=1, g=3, group elements intersect sink set partly.
        (7, [(0, 4), (1, 4), (2, 5), (3, 5), (4, 6), (5, 6)], (4, 5, 6),
         ((0, 1, 2), (3, 4, 5)), False),
        # 6. YES: d=3, t=1, g=3, G2 overlaps sinks.
        (9, [(0, 6), (1, 6), (2, 7), (3, 7), (4, 8), (5, 8)], (6, 7, 8),
         ((0, 1, 2), (3, 4, 5), (6, 7, 8)), True),
        # 7. YES: d=6, t=2, g=3 (3DM-like: two triples cover {6..11}, third conflicts).
        (12, [(0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11),
              (0, 9), (1, 10)], (6, 7, 8, 9, 10, 11),
         ((0, 1, 2), (3, 4, 5), (6, 7, 8)), True),
        # 8. NO: d=6, t=2, g=3 (all three triples share sink 6 or 7).
        (15, [(0, 9), (1, 10), (2, 11), (3, 9), (4, 12), (5, 13),
              (6, 10), (7, 13), (8, 14)], (9, 10, 11, 12, 13, 14),
         ((0, 1, 2), (3, 4, 5), (6, 7, 8)), False),
    ]
    # Add 6 seeded cyclic instances with d=3, t=1, g in {2, 3}.
    rng = Random(20261128)
    for _ in range(6):
        n0 = 9
        edges = [(u, v) for u in range(n0) for v in range(n0) if u != v and rng.random() < 0.25]
        sinks = (6, 7, 8)
        groups = ((0, 1, 2), (3, 4, 5))
        adj = [0] * n0
        for u, v in edges:
            adj[u] |= 1 << v
        yes_no = _solve_strict_gammoid_q_parity_exact(adj, sinks, groups) is not None
        test_cases.append((n0, edges, sinks, groups, yes_no))

    checked = 0
    for n0, edges, sinks, groups, expected_yes in test_cases:
        adj = [0] * n0
        for u, v in edges:
            adj[u] |= 1 << v
        parity_sol = _solve_strict_gammoid_q_parity_exact(adj, sinks, groups)
        assert (parity_sol is not None) == expected_yes

        rows, weights, k, vuln, threshold, meta = strict_gammoid_q_parity_to_actuator(
            adj, sinks, groups
        )
        n_star = len(rows)
        d = meta["d"]
        q = meta["q"]
        g_count = meta["g"]
        rho_star = matching_number(rows, (1 << n_star) - 1)
        assert n_star - rho_star == d
        assert k == q - 1

        # Verify source SCCs match the exact g group cycles.
        expected_cycles = frozenset(
            frozenset((s_i, *a_list, *f_list))
            for s_i, a_list, f_list in meta["group_vertices"]
        )
        actual_cycles = frozenset(
            frozenset(e for e in range(n_star) if mask >> e & 1)
            for mask in source_components(rows)
        )
        assert actual_cycles == expected_cycles

        # Verify loop and parallel copy identities in M*.
        for c_v in meta["centers"].values():
            assert dual_matching_rank(rows, (c_v,)) == 0
        for s_i, a_list, f_list in meta["group_vertices"]:
            assert dual_matching_rank(rows, (s_i,)) == 0
            for f_ir in f_list:
                assert dual_matching_rank(rows, (f_ir,)) == 0
        for i, (_, a_list, _) in enumerate(meta["group_vertices"]):
            for r_idx, a_ir in enumerate(a_list):
                v_ir = groups[i][r_idx]
                assert dual_matching_rank(rows, (a_ir,)) == dual_matching_rank(rows, (v_ir,))
                assert dual_matching_rank(rows, (a_ir, v_ir)) == dual_matching_rank(rows, (v_ir,))

        # Verify restricted dual matroid on {a_{i,r}} U Z matches M(D_1, T') on small cases.
        aff_nonloops = tuple(
            e for _, a_list, _ in meta["group_vertices"] for e in a_list
        ) + meta["z_vertices"]
        proj_map = {
            a_list[r_idx]: groups[i][r_idx]
            for i, (_, a_list, _) in enumerate(meta["group_vertices"])
            for r_idx in range(q)
        }
        for z_v in meta["z_vertices"]:
            proj_map[z_v] = z_v
        if len(aff_nonloops) <= 11:
            for sub in subsets(aff_nonloops):
                r_star = dual_matching_rank(rows, sub)
                r_d1 = max_vertex_disjoint_paths(
                    meta["d1_adj"], (proj_map[e] for e in sub), meta["t_prime"]
                )
                assert r_star == r_d1

        # Verify exact actuator threshold decision across all affordable group+Z choices.
        found_feasible = False
        for mask_g in range(1 << g_count):
            chosen = set()
            for i, (s_i, a_list, _) in enumerate(meta["group_vertices"]):
                if mask_g >> i & 1:
                    chosen.update(a_list)
                else:
                    chosen.add(s_i)
            for z_sub in subsets(meta["z_vertices"]):
                cand = frozenset(chosen | z_sub)
                cost = sum(weights[e] for e in cand)
                if cost > threshold:
                    continue
                witness = frozenset(cand & set(aff_nonloops))
                if verify_placement(rows, k, vuln, cand, witness):
                    # Verify full graph survival (reachability + matching) under all <= k failures.
                    full_star = (1 << n_star) - 1
                    for fail in subsets(cand & vuln, k):
                        surv = cand - fail
                        reached = set(surv)
                        front = list(surv)
                        for v in front:
                            row_m = rows[v]
                            while row_m:
                                bit = row_m & -row_m
                                row_m ^= bit
                                u = bit.bit_length() - 1
                                if u not in reached:
                                    reached.add(u)
                                    front.append(u)
                        assert len(reached) == n_star
                        surv_mask = sum(1 << e for e in surv)
                        assert matching_number(rows, full_star ^ surv_mask) == n_star - len(surv)
                    found_feasible = True
        assert found_feasible == expected_yes
        checked += 1

    return checked


if __name__ == "__main__":
    d012_count, d_counts, large_d1 = audit_small_deficiency()
    print(
        f"PASS: {d012_count} exhaustive d=0,1,2 optimum checks "
        f"(counts {d_counts}) and {large_d1} larger d=1 checks."
    )
    ip_cnt, mimic_cnt, strict_cnt, kk_info = audit_gammoid_and_cut_covering()
    print(
        f"PASS: {ip_cnt} Ingleton-Piff rank tables, {mimic_cnt} mimicking graphs "
        f"({strict_cnt} with rho_Z < rho), and {len(kk_info)} Koana-Kumabe F_p checks."
    )
    sp_cnt, ax_cnt, pair_cnt = audit_sparse_paving_lower_bound()
    print(
        f"PASS: {sp_cnt} sparse-paving oracle instances, {ax_cnt} matroid axiom checks, "
        f"and {pair_cnt} indistinguishable-pair optima checked."
    )
    qp_cnt = audit_strict_gammoid_q_parity_reduction()
    print(f"PASS: {qp_cnt} Strict Gammoid 3-Parity actuator reductions checked.")
