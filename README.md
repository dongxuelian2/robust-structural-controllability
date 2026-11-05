# Robust Structural Controllability under Actuator Failures

Mathematical research on fault-tolerant actuator placement in structured linear systems. The focus is on approximation algorithms, exact algorithms for restricted digraph classes, and verifiable combinatorial bounds.

## Research target

For the one-input-failure robust minimal structural controllability problem (1-RMSCP), study whether the known **5/2-approximation** for general directed graphs can be improved to **2**. The 2026 reference already provides a **2-approximation for strongly connected digraphs**.

Secondary questions:

- Exact algorithms or structural characterizations for restricted classes, including bounded-treewidth and planar digraphs.
- Tightness examples and improved lower bounds for existing approximation schemes.
- Extensions to multiple input failures when supported by new mathematical results.

Distinguish the dedicated-input and unrestricted-input models, including whether multiple inputs may actuate the same state. All claims should use the precise model and cost measure from the cited paper.

## Starting reference

Diogo Poças, Sérgio Pequito, and Guilherme Ramos, **On the NP-hardness of the robust minimal structural controllability problem**, *Automatica*, 190 (2026), article 113074.

DOI: https://doi.org/10.1016/j.automatica.2026.113074

The reference establishes NP-hardness even for one input failure and approximation guarantees of 5/2 in general digraphs and 2 in strongly connected digraphs.

## Research outputs

Keep proved theorems, precise counterexamples, and minimal exact computational checks in this repository. A self-contained theorem for a restricted graph class is useful even if the main approximation target is not attained.

[Degree threshold and sparse hardness](hardness-degree-threshold.md): exact optimization for state outdegree at most two; APX-hardness with degrees in {2,3} and no mandatory actuator, including strongly connected bipartite graphs with undirected girth at least N^(1−delta) for every fixed 0<delta<1. The reduction preserves OPT=2|E|+tau and gives a 2238/2233−epsilon inapproximability bound using the absolute published Vertex Cover gap. [Reduction obstructions](hardness-obstructions.md) explain the failure of one-leaf compression and a 5/4 barrier for forced-basis covering constructions.

[Source-SCC covering reduction](source-scc-cover-hardness.md): pairing private basis elements in source two-cycles changes the cost relation to `OPT=2*(|U|+tau)` and gives **250/249−epsilon** hardness for general directed graphs. Source-loop completion removes all mandatory locations; exact odd padding gives undirected girth at least N^(1−delta). The targets remain loopless, bipartite and subcubic, with maximum indegree and outdegree three. All added states remain selectable. Unlike the preceding class, these graphs have several source SCCs; the two-copy construction preserves both solution directions without parity exceptions.

[Variable-basis path-cover equivalence](variable-basis-path-cover.md): an unforced construction with exact optimum `n+pc(H)`, polynomial solution maps, degree-three compression and the same girth padding. The classical (1,2)-TSP 8/7 algorithm limits this particular route to stronger hardness. [Private-basis hyperforest objective](private-basis-hyperforest.md) extends the construction to larger neighborhoods, proves `OPT=n+min_C(|C|-r_C)`, and explains why frequency-two and regular incidence sources lose their covering gaps.

[Sparse cotransversal transfer](cotransversal-hardness-transfer.md): every coloop-free positive-rank cotransversal matroid with a full-rank presentation of row size at most three has a strongly connected subcubic control realization with exactly the same robust optimum, polynomial solution maps and the same girth padding. A certified example is nontransversal for every basis; an explicit strongly connected source family refutes a head/tail-gate reduction from directed path cover. An acyclic gammoid-completion family has the sharp rank-preserving distortion 2r/(r+1), even though every added resource is individually avoidable in an optimum.

[Strongly connected degree universality](row-splitting-equivalence.md) removes the row-size restriction: every strongly connected directed instance has a loopless, bipartite, bidirected realization with degrees in {2,3}, no mandatory actuator, exactly the same OPT, cost-nonincreasing polynomial solution maps and the same girth padding. The underlying exchange applies to every coloop-free positive-rank cotransversal matroid. [Uniform assignment transfer](assignment-hardness-transfer.md) additionally handles arbitrary unit-cost uniform robust assignments without an actuator baseline; a dominated-resource exchange shows that adding a private basis preserves OPT and extends the hyperforest formula to sources without private resources.

[Basis-relaxation obstruction](basis-relaxation-gap.md): replacing a matroid by the assignment model of its fundamental exchanges at an arbitrary basis can underestimate robust OPT by a factor approaching 3/2, even for simple, cosimple fundamental transversal matroids with resources serving at most two jobs. Maximizing total column degrees in a maximal dual presentation does not fix the loss. This is a proved relaxation gap, not an NP-hardness claim. All 513 bases of the nontransversal prototype are also audited; every basis contained in an optimum placement has a strict relaxation gap in that example.

Exact checks (Python 3.10+, standard library only): `python verify_outdegree_two.py`, `python verify_hardness_gadgets.py`, `python audit_hardness.py`, `python verify_variable_basis.py`, `python verify_private_basis_spares.py`, `python verify_gammoid_obstructions.py`, `python verify_assignment_transfer.py`, `python verify_row_splitting.py`, and `python verify_basis_relaxation.py`.
