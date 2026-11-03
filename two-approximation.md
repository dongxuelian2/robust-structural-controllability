# A quota-first factor-two approximation for general digraphs

## Scope and literature status

Model: one dedicated actuator per selected state, no duplicate placements, and robustness against deletion of any one actuator. The state graph is nonempty. All costs are one. The method below applies to arbitrary directed state graphs, including self-loops. It asserts infeasibility precisely when a source SCC is a singleton.

The starting reference is Poças, Pequito and Ramos, *Automatica* 190 (2026), 113074, DOI https://doi.org/10.1016/j.automatica.2026.113074. Its published abstract states a 5/2 approximation for general digraphs and a factor-two approximation for strongly connected digraphs. This note proves a factor-two guarantee for general digraphs under the dedicated, no-duplicates model. Its literature novelty has not been established. Related matroid background is Matthias Bentert, Fedor V. Fomin, Petr A. Golovach and Laure Morelle, *Fault-Tolerant Matroid Bases*, ESA 2025, DOI https://doi.org/10.4230/LIPIcs.ESA.2025.83. Robust spanning and bounded fault-tolerant bases should not be presented as newly introduced concepts here. The specific claim to assess against prior work is the quota-first relaxation and its general-digraph approximation guarantee.

## Notation

Let T be the transversal matroid of matchable subsets of the right-side state vertices; its rank is ν. Write M=T*, d=n−ν, and source SCCs C_1,...,C_ℓ. For any selected actuator set P, matching feasibility is equivalent to P spanning M, equivalently V\P independent in T. Ordinary structural controllability also requires P to hit each source SCC.

Single-failure robustness is equivalent to P−p spanning M for every p∈P, together with |P∩C_i|≥2 for every source SCC. A spanning set's fragile elements are exactly the coloops of the restricted matroid M|P.

## Algorithm

1. Compute source SCCs. If any is a singleton, return infeasible.
2. Compute a minimum-cardinality set P satisfying (a) P spans M, and (b) |P∩C_i|≥2 for every i. Let L=|P|. This is a relaxation of robust placement, not the ordinary minimal controllability problem.
3. Find the fragile elements D={b∈P: P−b does not span M}.
4. For each b∈D choose one e_b∉P such that P−b+e_b spans M. Return S=P∪{e_b:b∈D}.

All choices in step 4 can be made against the same original P. Duplicated repair states are installed only once.

## Exact solution of step 2 by a maximum flow

Write I=V\P. The requirements become:

- I is matchable by state edges;
- |I∩C_i|≤|C_i|−2 for each source SCC.

Create a network with a source s, one left and one right copy of each state, one group node per source SCC, and a sink t. Give s→left edges capacity one, state edges left→right capacity one, and each right vertex an outgoing edge of capacity one. For a right vertex in C_i, that edge leads to group i, whose edge to t has capacity |C_i|−2. For a vertex in no source SCC, it leads directly to t. An integral maximum flow defines a matching covering a maximum-cardinality I with the required quotas. Thus P=V\I is an optimum relaxation solution. This is also the intersection of the transversal matroid T with a partition matroid, but no general matroid-intersection implementation is needed.

## Repair existence and size

Every vertex has an incoming edge when every source SCC has size at least two. Hence T has no loops and M has no coloops: V−b spans M for every b.

For b∈D, P−b has rank d−1. Since V−b has rank d, some e_b∈V\P raises the rank to d; this proves repair existence.

The set D of coloops of M|P is independent in M, so |D|≤d. More explicitly, every basis B⊆P contains every element of D, and |B|=d.

The returned S is robust:

- If a repair-state actuator fails, the original spanning set P remains.
- If b∈D fails, the spanning set P−b+e_b remains.
- If b∈P\D fails, the spanning set P−b remains.

Source SCCs still have at least one actuator after any failure because P already had two. Thus S is a valid robust placement and |S|≤L+|D|≤L+d.

## Approximation guarantee and certificates

Let OPT be the minimum robust cardinality. Every robust solution is feasible for step 2, so L≤OPT. If d>0, removing any actuator from a robust solution must leave at least d states to span M, hence OPT≥d+1. Consequently

|S| ≤ L+d ≤ OPT+d ≤ 2 OPT−1 < 2 OPT.

For d=0, D is empty, and the algorithm is exact with L=2ℓ. This proves a factor-two approximation on all feasible digraphs.

The directly computable lower bound max(L,d+1) for d>0 (or L for d=0), together with the returned |S|, yields an a posteriori ratio certificate. Stronger instance-specific bound: |S|≤L+|D|.

## Cubic implementation without repeated matching computations

Retain the matching F from step 2, covering I=V\P. For each b∈P, look for an augmenting path starting at the currently unmatched right vertex b, traversing an incident state edge backward to a left vertex and then its matching edge forward to the matched right vertex. If an unmatched left vertex is reached, I+b is matchable and b is not fragile. Otherwise b∈D.

For a fragile b, choose any incoming state edge u→b. Its left endpoint u must be matched by F; otherwise there would be an augmenting path of length one. Let u→e_b be its matching edge. Because e_b∈I, swapping u→e_b for u→b gives a matching covering I−e_b+b. Therefore P−b+e_b spans M. This gives the repair directly, without a new matching call.

There are O(n) network vertices and O(n²) edges, integral maximum-flow value at most n, and each augmenting search takes O(n²). A straightforward augmenting-path flow implementation therefore takes O(n³). Each of at most n fragility searches also takes O(n²). SCC computation fits within this bound. Total running time O(n³), space O(n²), under the standard graph RAM model. The reference bitmask representation does not alter the combinatorial operation bound, though arbitrary-precision bit operations have their usual bit-complexity costs.

## Files and verification

`two_approx_solver.py` is the self-contained algorithm, apart from reusing `source_components` from `bounded_core_solver.py`. Its output includes the relaxation placement, relaxation value, critical states, repair states, and a matching certificate.

`check_two_approx.py` brute-forces all actuator sets on each tested graph, independently verifying exact relaxation optimality, robust feasibility of the returned witness, |D|≤d, and the approximation inequalities. The stated proof is independent of finite testing. Running `python check_two_approx.py` passed all 66,066 directed graphs with self-loops on one through four vertices and 600 fixed-seed random graphs on five through eight vertices: 66,666 total, 38,951 feasible. The largest observed ratio was 1.5, which is only an experimental observation and not a proved approximation factor. Run `python audit_cubic_solver.py` for the independent implementation and oracle audit, which separately passed 66,066 exhaustive and 400 random instances.

## Limitations

The proof does not establish an approximation ratio strictly bounded below two uniformly in n, optimality of this approximation factor, weighted-cost guarantees, multiple-failure guarantees, or publication novelty. In particular, |S|≤L+d is an additive statement in unit costs and cannot simply be applied to arbitrary actuator costs.

## Asymptotic tightness of the independent-repair strategy

For each d≥2, take vertices c_i,e_i,f_i (i=0,...,d−1), h, and g. Add edges c_i→e_i,f_i; e_i,f_i→c_i,h; h→e_i,f_i,g for every i; and g→c_i,h for every i. The graph is strongly connected and has n=3d+2 states.

Only the d vertices c_i and h have outgoing edges into the 2d+1 target vertices {e_i,f_i,g}, so at most d+1 targets are matched. At most d+1 center vertices {c_i,h} are matched. The matching F={c_i→e_i, e_i→c_i: all i}∪{h→f_0,f_0→h} attains 2d+2 edges. Thus the deficiency is exactly d.

The set W={e_0,...,e_{d−1},g} is robust: after deleting g, match f_i via c_i for all i and g via h; after deleting e_j, match f_i via c_i and the additional unactuated target e_j via h. In both cases match every center c_i via the left vertex e_i and h via the left vertex g. These left endpoints are distinct from c_i and h, so all these edges coexist in a matching. Strong connectivity gives accessibility after every deletion. Hence OPT≤d+1; the rank lower bound proves OPT=d+1.

The maximum-flow relaxation may return F above, whose unmatched set is P={f_1,...,f_{d−1},g}, with |P|=d; this is optimal for the relaxation. It is a basis of M, so every one of its d elements is fragile. Choosing incoming c_i→f_i repairs f_i with e_i for i≥1; choosing h→g repairs g with f_0. These repairs are distinct, so the returned set has 2d elements, and its approximation ratio is 2d/(d+1), tending to two.

Thus no uniform guarantee strictly below two holds for the strategy with arbitrary optimum-relaxation and independent-repair choices. This is not an inapproximability result for the original problem, and does not preclude better tie-breaking, joint repair, or postprocessing. `check_two_approx_tightness.py` verifies the displayed optimal witness and the actual reference implementation's 2d output for d=2,...,25.
