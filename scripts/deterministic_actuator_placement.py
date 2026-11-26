"""Deterministic algorithms and exact barriers for weighted actuator placement."""

from itertools import combinations
from math import inf

from bounded_core_solver import matching_number, source_components
from linear_accessibility_quotas import PurchaseSolution, solve_represented
from weighted_fault_core import weighted_input


def maximum_matching(rows, allowed):
    """Return a maximum matching size and right-to-left partner map for B_G."""
    n = len(rows)
    mate = [-1] * n
    count = 0
    for start in range(n):
        queue = [start]
        seen_left = {start}
        predecessor = {}
        free = None
        for left in queue:
            mask = rows[left] & allowed
            while mask:
                bit = mask & -mask
                mask ^= bit
                right = bit.bit_length() - 1
                if right in predecessor:
                    continue
                predecessor[right] = left
                if mate[right] == -1:
                    free = right
                    break
                other = mate[right]
                if other not in seen_left:
                    seen_left.add(other)
                    queue.append(other)
            if free is not None:
                break
        if free is None:
            continue
        left_mate = {left: right for right, left in enumerate(mate) if left >= 0}
        right = free
        while True:
            left = predecessor[right]
            old = left_mate.get(left)
            mate[right] = left
            if old is None:
                break
            right = old
        count += 1
    return count, tuple(mate)


def _validate_state_graph(rows, weights, k, vulnerable):
    weights, ground, vulnerable = weighted_input(weights, k, vulnerable)
    rows = tuple(rows)
    n = len(rows)
    if n != len(weights) or n == 0 or any(
        type(row) is not int or not 0 <= row < 1 << n for row in rows
    ):
        raise ValueError("Expected a nonempty graph and one weight per state.")
    return rows, weights, ground, vulnerable


def strict_gammoid_digraph(rows):
    """Build the Ingleton-Piff strict-gammoid digraph (adj, sinks, mate) for T^*."""
    rows = tuple(rows)
    n = len(rows)
    if n == 0 or any(type(row) is not int or not 0 <= row < 1 << n for row in rows):
        raise ValueError("Expected a nonempty graph of n-bit out-neighbor masks.")
    _, mate = maximum_matching(rows, (1 << n) - 1)
    sinks = tuple(r for r, left in enumerate(mate) if left == -1)
    adj = [0] * n
    for r, left in enumerate(mate):
        if left != -1:
            adj[r] = rows[left] & ~(1 << r)
    return tuple(adj), sinks, mate


def solve_deficiency_zero(rows, weights, k, vulnerable=None):
    """Solve weighted actuator placement exactly in polynomial time when d = 0."""
    rows, weights, ground, vulnerable = _validate_state_graph(
        rows, weights, k, vulnerable
    )
    n = len(rows)
    if matching_number(rows, (1 << n) - 1) != n:
        raise ValueError("solve_deficiency_zero requires matching deficiency d = 0.")
    k_eff = min(k, len(vulnerable))
    threshold = k_eff + 1
    selected = set()
    for mask in source_components(rows):
        group = frozenset(e for e in ground if mask >> e & 1)
        safe = sorted(group - vulnerable, key=lambda e: (weights[e], e))
        exposed = sorted(group & vulnerable, key=lambda e: (weights[e], e))
        options = [safe[:1]] if safe else []
        if len(exposed) >= threshold:
            options.append(exposed[:threshold])
        if not options:
            return None
        selected.update(
            min(options, key=lambda seq: (sum(weights[e] for e in seq), seq))
        )
    return PurchaseSolution(frozenset(selected), frozenset())


def solve_deficiency_one(rows, weights, k, vulnerable=None):
    """Solve weighted actuator placement exactly in polynomial time when d = 1."""
    rows, weights, ground, vulnerable = _validate_state_graph(
        rows, weights, k, vulnerable
    )
    n = len(rows)
    adj, sinks, _ = strict_gammoid_digraph(rows)
    if len(sinks) != 1:
        raise ValueError("solve_deficiency_one requires matching deficiency d = 1.")
    sink = sinks[0]

    # Nonloops of M = T^* are vertices that can reach the unique sink in D_H.
    rev = [0] * n
    for u, mask in enumerate(adj):
        while mask:
            bit = mask & -mask
            mask ^= bit
            v = bit.bit_length() - 1
            rev[v] |= 1 << u
    nonloop_mask = 1 << sink
    queue = [sink]
    for v in queue:
        pred_mask = rev[v] & ~nonloop_mask
        while pred_mask:
            bit = pred_mask & -pred_mask
            pred_mask ^= bit
            u = bit.bit_length() - 1
            nonloop_mask |= bit
            queue.append(u)

    k_eff = min(k, len(vulnerable))
    threshold = k_eff + 1
    safe_state = threshold + 1

    groups = tuple(
        frozenset(e for e in ground if mask >> e & 1)
        for mask in source_components(rows)
    )
    grouped = frozenset().union(*groups)
    ungrouped = ground - grouped

    # DP state: 0..threshold count vulnerable nonloops without a safe nonloop;
    # safe_state records selection of at least one safe nonloop.
    dp = {0: (0, ())}

    def relax(nxt, state, cost, chosen):
        prev = nxt.get(state)
        if prev is None or cost < prev[0] or (cost == prev[0] and chosen < prev[1]):
            nxt[state] = (cost, chosen)

    for group in groups:
        safe_nonloops = sorted(
            (e for e in group - vulnerable if nonloop_mask >> e & 1),
            key=lambda e: (weights[e], e),
        )
        safe_loops = sorted(
            (e for e in group - vulnerable if not nonloop_mask >> e & 1),
            key=lambda e: (weights[e], e),
        )
        vuln_nonloops = sorted(
            (e for e in group & vulnerable if nonloop_mask >> e & 1),
            key=lambda e: (weights[e], e),
        )
        vuln_loops = sorted(
            (e for e in group & vulnerable if not nonloop_mask >> e & 1),
            key=lambda e: (weights[e], e),
        )
        vn_pref = [0]
        for e in vuln_nonloops[:threshold]:
            vn_pref.append(vn_pref[-1] + weights[e])
        vl_pref = [0]
        for e in vuln_loops[:threshold]:
            vl_pref.append(vl_pref[-1] + weights[e])

        transitions = []
        if safe_nonloops:
            e = safe_nonloops[0]
            transitions.append((True, 0, weights[e], (e,)))
        if safe_loops:
            s = safe_loops[0]
            ws = weights[s]
            for a in range(min(threshold, len(vuln_nonloops)) + 1):
                transitions.append(
                    (False, a, ws + vn_pref[a], (s, *vuln_nonloops[:a]))
                )
        for a in range(min(threshold, len(vuln_nonloops)) + 1):
            need_vl = threshold - a
            if need_vl <= len(vuln_loops):
                transitions.append(
                    (
                        False,
                        a,
                        vn_pref[a] + vl_pref[need_vl],
                        (*vuln_nonloops[:a], *vuln_loops[:need_vl]),
                    )
                )
        if not transitions:
            return None

        nxt = {}
        for state, (base_cost, base_chosen) in dp.items():
            for to_safe, add_vn, delta, add_elems in transitions:
                if to_safe or state == safe_state:
                    nstate = safe_state
                else:
                    nstate = min(threshold, state + add_vn)
                relax(nxt, nstate, base_cost + delta, base_chosen + add_elems)
        dp = nxt

    safe_nonloops = sorted(
        (e for e in ungrouped - vulnerable if nonloop_mask >> e & 1),
        key=lambda e: (weights[e], e),
    )
    vuln_nonloops = sorted(
        (e for e in ungrouped & vulnerable if nonloop_mask >> e & 1),
        key=lambda e: (weights[e], e),
    )
    vn_pref = [0]
    for e in vuln_nonloops[:threshold]:
        vn_pref.append(vn_pref[-1] + weights[e])
    ungrouped_trans = [(False, 0, 0, ())]
    if safe_nonloops:
        e = safe_nonloops[0]
        ungrouped_trans.append((True, 0, weights[e], (e,)))
    for a in range(1, min(threshold, len(vuln_nonloops)) + 1):
        ungrouped_trans.append((False, a, vn_pref[a], tuple(vuln_nonloops[:a])))

    nxt = {}
    for state, (base_cost, base_chosen) in dp.items():
        for to_safe, add_vn, delta, add_elems in ungrouped_trans:
            if to_safe or state == safe_state:
                nstate = safe_state
            else:
                nstate = min(threshold, state + add_vn)
            relax(nxt, nstate, base_cost + delta, base_chosen + add_elems)
    dp = nxt

    candidates = [dp[s] for s in (threshold, safe_state) if s in dp]
    if not candidates:
        return None
    _, best_chosen = min(candidates, key=lambda pair: (pair[0], len(pair[1]), pair[1]))
    placement = frozenset(best_chosen)
    safe_witness = sorted(e for e in placement - vulnerable if nonloop_mask >> e & 1)
    if safe_witness:
        witness = frozenset((safe_witness[0],))
    else:
        vuln_witness = sorted(e for e in placement & vulnerable if nonloop_mask >> e & 1)
        witness = frozenset(vuln_witness[:threshold])
    return PurchaseSolution(placement, witness)


def rank_two_representation(rows):
    """Build an exact Q-representation with O(log n)-bit entries when d = 2."""
    rows = tuple(rows)
    n = len(rows)
    if n == 0 or any(type(row) is not int or not 0 <= row < 1 << n for row in rows):
        raise ValueError("Expected a nonempty graph of n-bit out-neighbor masks.")
    full = (1 << n) - 1
    rho = matching_number(rows, full)
    if n - rho != 2:
        raise ValueError("rank_two_representation requires matching deficiency d = 2.")
    nonloops = [
        e for e in range(n) if matching_number(rows, full ^ (1 << e)) == rho
    ]
    reps = []
    class_of = {}
    for e in nonloops:
        assigned = None
        for idx, rep in enumerate(reps):
            if matching_number(rows, full ^ (1 << e) ^ (1 << rep)) == rho - 1:
                assigned = idx
                break
        if assigned is None:
            assigned = len(reps)
            reps.append(e)
        class_of[e] = assigned
    columns = []
    for e in range(n):
        if e in class_of:
            columns.append((1, class_of[e] + 1))
        else:
            columns.append((0, 0))
    return tuple(columns)


def solve_deficiency_two(rows, weights, k, vulnerable=None):
    """Solve weighted actuator placement deterministically in FPT(k) when d = 2."""
    rows, weights, ground, vulnerable = _validate_state_graph(
        rows, weights, k, vulnerable
    )
    columns = rank_two_representation(rows)
    groups = tuple(
        frozenset(e for e in ground if mask >> e & 1)
        for mask in source_components(rows)
    )
    return solve_represented(columns, weights, k, vulnerable, groups)


def solve_small_deficiency(rows, weights, k, vulnerable=None):
    """Dispatch exact deterministic solver for matching deficiency d in {0, 1, 2}."""
    rows, weights, ground, vulnerable = _validate_state_graph(
        rows, weights, k, vulnerable
    )
    n = len(rows)
    rho = matching_number(rows, (1 << n) - 1)
    d = n - rho
    if d == 0:
        return solve_deficiency_zero(rows, weights, k, vulnerable)
    if d == 1:
        return solve_deficiency_one(rows, weights, k, vulnerable)
    if d == 2:
        return solve_deficiency_two(rows, weights, k, vulnerable)
    raise ValueError(f"Unsupported matching deficiency d = {d}; expected d <= 2.")


def mimicking_digraph(adj, sinks, terminal_subset, cut_cover):
    """Build torso digraph D_Z and dual bipartite graph B_Z from a cut cover Z.

    Vertices of augmented digraph D~_H are indexed as:
      - 0 .. n-1: original vertices V of adj
      - n .. n + |W| - 1: fresh source copies S = {w_in : w in W}
      - n + |W| .. n + |W| + d - 1: fresh sink copies T = {t_out : t in sinks}
    cut_cover is an iterable of vertex indices in D~_H containing S union T.
    """
    n = len(adj)
    terminal_list = sorted(terminal_subset)
    sinks = tuple(sinks)
    d = len(sinks)
    s_indices = tuple(range(n, n + len(terminal_list)))
    t_indices = tuple(range(n + len(terminal_list), n + len(terminal_list) + d))
    total = n + len(terminal_list) + d

    aug_adj = [set() for _ in range(total)]
    for u, mask in enumerate(adj):
        while mask:
            bit = mask & -mask
            mask ^= bit
            aug_adj[u].add(bit.bit_length() - 1)
    for pos, w in enumerate(terminal_list):
        aug_adj[n + pos].add(w)
    for pos, t in enumerate(sinks):
        aug_adj[t].add(n + len(terminal_list) + pos)

    z_set = frozenset(cut_cover)
    required = frozenset(s_indices) | frozenset(t_indices)
    if not required <= z_set or not z_set <= frozenset(range(total)):
        raise ValueError("cut_cover must lie in V(D~_H) and contain S union T.")

    # Order Z with Z \ T first so B_Z rows correspond to 0 .. |Z| - d - 1.
    t_set = frozenset(t_indices)
    z_non_t = sorted(z_set - t_set)
    z_t = sorted(t_set)
    z_ordered = tuple(z_non_t + z_t)
    z_pos = {v: i for i, v in enumerate(z_ordered)}

    torso_adj = [0] * len(z_ordered)
    for x_pos, x in enumerate(z_ordered):
        if x in t_set:
            continue
        seen = {x}
        queue = []
        for nxt in aug_adj[x]:
            if nxt in z_set:
                if nxt != x:
                    torso_adj[x_pos] |= 1 << z_pos[nxt]
            elif nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
        for u in queue:
            for nxt in aug_adj[u]:
                if nxt in z_set:
                    if nxt != x:
                        torso_adj[x_pos] |= 1 << z_pos[nxt]
                elif nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)

    # B_Z has left vertices Z \ T (matching M_0: z <-> z) and right vertices Z.
    b_z_rows = tuple(
        (1 << i) | torso_adj[i] for i in range(len(z_non_t))
    )
    index_map = {
        "z_ordered": z_ordered,
        "z_pos": z_pos,
        "s_in_b_z": tuple(z_pos[s] for s in s_indices),
        "t_in_b_z": tuple(z_pos[t] for t in t_indices),
        "w_to_b_z": {w: z_pos[n + pos] for pos, w in enumerate(terminal_list)},
    }
    return tuple(torso_adj), b_z_rows, index_map


def sparse_paving_lower_bound_instance(
    u, r, k, hidden_group_tuple=None, loopless=False
):
    """Build the sparse-paving rank-oracle lower bound instance M_H.

    When hidden_group_tuple is None, H contains all t-element subsets I of
    {0, ..., u - 1}; when given, H removes the single circuit-hyperplane H_{I*}.
    """
    if type(u) is not int or type(r) is not int or type(k) is not int or k < 1:
        raise ValueError("Expected integer parameters with k >= 1.")
    if loopless:
        if r < 3:
            raise ValueError("Loopless variant requires r >= 3.")
        t = (r - 1) // 2
        has_parity_coloop = (r - 1) % 2 == 1
    else:
        if r < 2:
            raise ValueError("Primary variant requires r >= 2.")
        t = r // 2
        has_parity_coloop = r % 2 == 1
    if u < t + 1:
        raise ValueError(f"Expected u >= t + 1 = {t + 1}.")

    m_cost = 2 * t + 5
    weights = []
    vulnerable = set()
    groups = []
    pair_ids = []
    loop_like_ids = []

    for _ in range(u):
        g = set()
        s_id = len(weights)
        weights.append((k + 1) * m_cost)
        g.add(s_id)
        loop_like_ids.append(s_id)

        a_id = len(weights)
        weights.append(m_cost + 1)
        vulnerable.add(a_id)
        g.add(a_id)

        b_id = len(weights)
        weights.append(m_cost + 1)
        vulnerable.add(b_id)
        g.add(b_id)
        pair_ids.append((a_id, b_id))

        for _ in range(k - 1):
            f_id = len(weights)
            weights.append(m_cost)
            vulnerable.add(f_id)
            g.add(f_id)
            loop_like_ids.append(f_id)

        groups.append(frozenset(g))

    z_ids = []
    for _ in range(k):
        z_id = len(weights)
        weights.append(0)
        vulnerable.add(z_id)
        z_ids.append(z_id)

    c0_id = None
    if has_parity_coloop:
        c0_id = len(weights)
        weights.append(0)

    clp_id = None
    if loopless:
        clp_id = len(weights)
        weights.append(0)
        loop_like_ids.append(clp_id)

    core_ids = frozenset(e for pair in pair_ids for e in pair) | frozenset(z_ids)
    lp_class = frozenset(loop_like_ids)

    hidden_set = None
    if hidden_group_tuple is not None:
        hidden_tuple = tuple(sorted(hidden_group_tuple))
        if len(hidden_tuple) != t or len(set(hidden_tuple)) != t or any(
            not 0 <= i < u for i in hidden_tuple
        ):
            raise ValueError("Invalid hidden_group_tuple.")
        hidden_set = frozenset(
            e for i in hidden_tuple for e in pair_ids[i]
        )

    # Precompute pair mask lookup: a subset of S_core is in H_all iff it has size 2t
    # and is a union of t full {a_i, b_i} pairs.
    elem_to_pair_idx = {}
    for idx, (a_id, b_id) in enumerate(pair_ids):
        elem_to_pair_idx[a_id] = idx
        elem_to_pair_idx[b_id] = idx

    def is_circuit_hyperplane(core_sub):
        if len(core_sub) != 2 * t:
            return False
        if hidden_set is not None and core_sub == hidden_set:
            return False
        counts = {}
        for e in core_sub:
            idx = elem_to_pair_idx.get(e)
            if idx is None:
                return False
            counts[idx] = counts.get(idx, 0) + 1
        return len(counts) == t and all(c == 2 for c in counts.values())

    def rank_fn(subset):
        subset = frozenset(subset)
        core_sub = subset & core_ids
        r_core = 2 * t - 1 if is_circuit_hyperplane(core_sub) else min(len(core_sub), 2 * t)
        extra = 0
        if c0_id is not None and c0_id in subset:
            extra += 1
        if loopless and (subset & lp_class):
            extra += 1
        return r_core + extra

    opt_yes = (k + 1) * m_cost * u + 2 * t
    opt_no = (k + 1) * m_cost * u + 2 * t + 2
    return (
        len(weights),
        rank_fn,
        tuple(weights),
        frozenset(vulnerable),
        tuple(groups),
        opt_yes,
        opt_no,
        hidden_set,
    )


def strict_gammoid_q_parity_to_actuator(adj, sinks, groups):
    """Reduce Strict Gammoid q-Parity to weighted actuator placement preserving d = q*t."""
    adj = tuple(adj)
    n0 = len(adj)
    sinks = tuple(sinks)
    d = len(sinks)
    groups = tuple(tuple(g) for g in groups)
    if not groups:
        raise ValueError("Expected at least one group.")
    q = len(groups[0])
    if q < 3 or any(len(g) != q for g in groups):
        raise ValueError("Expected all groups to have uniform size q >= 3.")
    if d == 0 or d % q != 0:
        raise ValueError("Expected |T| = d to be a positive multiple of q.")
    seen_grouped = set()
    for g in groups:
        if any(not 0 <= v < n0 for v in g) or seen_grouped.intersection(g):
            raise ValueError("Groups must be pairwise disjoint subsets of V(D).")
        seen_grouped.update(g)
    if len(set(sinks)) != d or any(not 0 <= t_v < n0 for t_v in sinks):
        raise ValueError("Sinks must be distinct vertices of V(D).")

    g_count = len(groups)
    m_big = 2 * d + 1
    forbidden_weight = 2 * (q * m_big + d + 1)
    k = q - 1

    # D_1 vertices: V (0..n0-1), T' (n0..n0+d-1), Z (n0+d..n0+d+q-2).
    n_tprime = n0
    n_z = n0 + d
    n1 = n0 + d + (q - 1)
    t_prime_set = set(range(n_tprime, n_tprime + d))
    v1_non_tprime = [v for v in range(n1) if v not in t_prime_set]

    d1_adj = [0] * n1
    for u in range(n0):
        d1_adj[u] = adj[u] & ~(1 << u)
    for pos, v in enumerate(sinks):
        d1_adj[v] |= 1 << (n_tprime + pos)
    t_prime_mask = sum(1 << v for v in t_prime_set)
    for z in range(n_z, n_z + q - 1):
        d1_adj[z] = t_prime_mask

    # Assign state indices in G*:
    # - V_1: 0 .. n1-1
    # - c_v for v in V_1 \ T': next |V_1 \ T'| vertices
    # - Group i (0 .. g_count-1): s_i, a_{i,1..q}, f_{i,0..q}
    center_of = {}
    cursor = n1
    for v in v1_non_tprime:
        center_of[v] = cursor
        cursor += 1

    group_vertices = []
    for _ in range(g_count):
        s_i = cursor
        cursor += 1
        a_list = tuple(range(cursor, cursor + q))
        cursor += q
        f_list = tuple(range(cursor, cursor + q + 1))
        cursor += q + 1
        group_vertices.append((s_i, a_list, f_list))

    total = cursor
    rows = [0] * total
    weights = [forbidden_weight] * total
    vulnerable = set(range(total))

    for z in range(n_z, n_z + q - 1):
        weights[z] = 0

    for v in v1_non_tprime:
        c_v = center_of[v]
        rows[v] = 1 << c_v
        rows[c_v] = (1 << v) | d1_adj[v]

    for i, (s_i, a_list, f_list) in enumerate(group_vertices):
        weights[s_i] = q * m_big
        vulnerable.remove(s_i)
        for a_ir in a_list:
            weights[a_ir] = m_big + 1

        rows[s_i] = 1 << f_list[0]
        for r in range(q):
            rows[f_list[r]] = (1 << a_list[r]) | (1 << groups[i][r])
            rows[a_list[r]] = 1 << f_list[r + 1]
        rows[f_list[q]] = 1 << s_i

    # Loop-only connector from f_{1, q} to every center c_v.
    f_1_q = group_vertices[0][2][q]
    for v in v1_non_tprime:
        rows[f_1_q] |= 1 << center_of[v]

    threshold = q * m_big * g_count + d
    metadata = {
        "n1": n1,
        "t_prime": tuple(range(n_tprime, n_tprime + d)),
        "z_vertices": tuple(range(n_z, n_z + q - 1)),
        "centers": center_of,
        "group_vertices": tuple(group_vertices),
        "d1_adj": tuple(d1_adj),
        "d": d,
        "q": q,
        "g": g_count,
        "m_big": m_big,
    }
    return tuple(rows), tuple(weights), k, frozenset(vulnerable), threshold, metadata
