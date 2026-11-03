# Exact bounded-core reduction for arbitrary matching deficiency

## Model and status

A directed state graph G=(V,E), n>0, one dedicated actuator at each selected state S⊆V, no repeated actuators at one state. Any one selected actuator may fail. Edges and remaining input parameters are independently free structural parameters. This note proves an exact reduction and a polynomial-time algorithm for every fixed matching deficiency d. It does not claim a new result in the literature, an FPT algorithm, or polynomial complexity when d is unbounded. Publication novelty requires literature review.

Let T be the transversal matroid on right-side state vertices of the state bipartite graph: independent subsets are those covered by some matching. Its rank is ν, the maximum matching size. Let M=T* be its dual and d=n−ν. For X⊆V,

r_M(X)=|X|+ν(V\X)−ν(V),

where ν(Y) is the largest matching restricted to right vertices Y. Let C_1,...,C_ℓ be the source strongly connected components.

A set S is structurally controllable iff r_M(S)=d and S meets each C_i. Indeed a dedicated-input matching covers all states iff every unactuated state can be matched using state edges, i.e. V\S is independent in T, equivalent to S spanning M. Accessibility is equivalent to meeting every source SCC.

Consequently S is robust to one actuator failure iff

1. r_M(S\{s})=d for every s∈S, and r_M(S)=d;
2. |S∩C_i|≥2 for every i.

The explicit no-failure condition prevents vacuous acceptance of the empty set.

## Feasibility

A robust solution exists iff every source SCC has at least two vertices.

Necessity follows from condition 2. For sufficiency, these SCC sizes imply every state has an incoming edge: vertices in a nontrivial source SCC have an internal incoming edge; vertices in non-source SCCs either have one internally or, for singleton SCCs, externally. Thus T has no loops, so M has no coloops. Selecting all states therefore remains spanning in M after any single removal, and also satisfies accessibility. This also follows directly: after removing the actuator at v, an incoming state edge u→v plus dedicated actuators at all remaining vertices gives a state-covering matching, and source accessibility remains.

## Bounded robust-core lemma

For any rank-d matroid M, every spanning set S whose removal of any one element remains spanning contains a set K⊆S with |K|≤2d having the same property.

Proof. Choose a basis B⊆S, |B|=d. For each b∈B, because S\{b} spans, the rank-(d−1) set B\{b} can be augmented by an element e_b∈S\B. Thus B−b+e_b is a basis. Put K=B∪{e_b:b∈B}. There are at most 2d elements. Removing b∈B leaves the indicated replacement basis. Removing any e∈K\B leaves B. Thus every one-element removal spans. For d=0 use B=K=∅. ∎

## Exact optimization formula

Let R_d={K⊆V: |K|≤2d, r_M(K)=d, and r_M(K−k)=d for every k∈K}. If some source SCC is a singleton, OPT=∞. Otherwise

OPT = min_{K∈R_d} ( |K| + Σ_i max(0, 2−|K∩C_i|) ).

Proof, upper bound. Complete a candidate K by adding exactly max(0,2−|K∩C_i|) unused vertices from each source SCC. This is possible because |C_i|≥2, and source SCCs are disjoint. Robust spanning is upward closed: if a newly added vertex fails, K survives; if a core vertex fails, its remaining spanning set survives. Thus the completion is robust and has the stated size.

Lower bound. Any optimal solution S contains a robust core K of size at most 2d by the lemma. To meet each source SCC twice, S\K contains at least max(0,2−|K∩C_i|) distinct vertices there. Sum over disjoint source SCCs. This yields the stated objective at K no greater than |S|. ∎

Algorithm: compute source SCCs and ν, enumerate subsets of sizes ≤min(n,2d), evaluate their rank and single-deletion ranks through matching calls, and minimize the completion cost. Number of candidates Σ_{j=0}^{min(n,2d)} binom(n,j). At most 2d+1 rank-oracle calls per candidate. This is an XP algorithm parameterized by d, polynomial for fixed d, not an FPT guarantee. It returns both an optimum count and an explicit actuator set.

## Recovery of d≤1 formulas

For d=0 the sole needed core is ∅, hence OPT=2ℓ on feasible instances.

For d=1, let U be the nonloops of M, equivalently states exposed by some maximum state matching of size n−1. The minimal robust cores are precisely pairs of distinct elements of U. Put q=Σ_i min(2,|U∩C_i|). Such a core can place min(2,q) of its vertices in source SCCs, so OPT=2ℓ+max(0,2−q). Feasibility guarantees at least two nonloops in M.

## Deficiency two: explicit classification

In rank two, discard loops and partition nonloops into parallel classes: x,y lie in the same class iff r_M({x,y})=1. A spanning restriction has no coloops precisely when it either meets at least three parallel classes or meets exactly two classes with at least two elements in each. Therefore every inclusion-minimal robust core has one of the forms:

- Three vertices in three distinct parallel classes.
- Four vertices, two in each of two distinct parallel classes.

Conversely every such triple or quadruple is robust spanning. Substitute these cores in the exact formula. Rank preprocessing needs singleton and pair queries; candidate enumeration then takes O(n^4) elementary candidate evaluations, or O(n^4 ℓ) with a naive source-count implementation. No higher-order matching checks are needed after class preprocessing. Source contributions can be evaluated in constant time per candidate after assigning each vertex its source-SCC label (or none).

Proof of classification. Rank-two spanning requires two classes. If exactly two occur, removing a unique member of either class reduces rank, so both multiplicities must be at least two; two from each suffice. If at least three occur, one from any three suffices. Minimality excludes any extra elements or loops. ∎

## Tightness and a counterexample to a naive closed form

For any d≥1 construct vertices a_i,b_i,c_i, i=1,...,d. Add c_i→a_i, c_i→b_i, a_i→c_i, b_i→c_i. For d≥2 connect centers in a directed cycle c_i→c_{i+1}, indices modulo d (for d=1 no extra center edge is needed). The graph is strongly connected, ℓ=1. At most d leaf vertices and d centers can be matched, and 2d is achievable, so its matching deficiency is d. Every maximum matching matches all centers and exactly one leaf from each pair {a_i,b_i}. Thus M is the direct sum of d rank-one two-element classes, with the center vertices loops. Robust spanning requires selecting both leaves from every pair: OPT=2d. This proves the core bound 2d is tight even for strongly connected state graphs.

For d=2 this gives a six-vertex strongly connected graph with OPT=4. By contrast the bidirected star with one center and three leaves has d=2, ℓ=1, dual matroid U_{2,3} on its leaves, and OPT=3. Thus d and ℓ alone do not determine OPT; in particular max(d+1,2ℓ) is not always the optimum. The parallel-class information genuinely matters beyond d=1.

## Validation

Use Python 3.10 or newer. Run `python verify_core_theorem.py`. It checks all directed graphs with loops on n=1,2,3 (530 instances), all loopless directed graphs on n=4 (4,096 instances), and 500 fixed-seed random graphs on n=4,...,8. Every actuator subset is brute-forced for each graph. Checks cover feasibility, bounded-core formula, d=0 and d=1 formulas, and the rank-two minimal-core classification. For the first 530 graphs an independent augmented-matching dynamic-programming oracle checks structural controllability before and after each actuator deletion. Finite testing supplements, and does not replace, the proofs above.

Reference solver: `bounded_core_solver.py` implements the general exact theorem with iterative (nonrecursive) matching, returns an explicit optimum placement, rejects singleton source SCCs, and validates input bitmasks. Run `python audit_bounded_solver.py` for an independent subset-DP matching oracle and per-vertex DFS SCC oracle. This audit checks all 66,066 directed graphs with loops on n=1,...,4, plus 400 fixed-seed random graphs on n=5,...,8 (seed 20261103). It verifies the optimum, infeasibility, and every returned placement under each actuator deletion, as well as input nonmutation and invalid-input rejection. The solver uses no external packages.

Recorded main validation: 5,126 graphs passed, deficiency histogram {0:1904, 1:2398, 2:667, 3:99, 4:40, 5:16, 6:2}. Explicit tight-family cases d=1,2,3 returned optima 2,4,6, respectively; the four-state bidirected star returned deficiency 2 and optimum 3.
