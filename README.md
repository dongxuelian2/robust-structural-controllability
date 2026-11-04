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

[Degree threshold and sparse hardness](hardness-degree-threshold.md): exact optimization for state outdegree at most two; APX-hardness at degree three, including strongly connected bipartite graphs with undirected girth at least N^(1−delta) for every fixed 0<delta<1. The reduction preserves OPT=2|E|+tau and gives a 469/468−epsilon inapproximability bound. [Reduction obstructions](hardness-obstructions.md) explain the failure of one-leaf compression and a 5/4 barrier for forced-basis covering constructions.

[Variable-basis path-cover equivalence](variable-basis-path-cover.md): an unforced construction with exact optimum `n+pc(H)`, polynomial solution maps, degree-three compression and the same girth padding. The classical (1,2)-TSP 8/7 algorithm limits this particular route to stronger hardness. [Private-basis hyperforest objective](private-basis-hyperforest.md) extends the construction to larger neighborhoods, proves `OPT=n+min_C(|C|-r_C)`, and explains why frequency-two and regular incidence sources lose their covering gaps.

Exact checks (Python 3.10+, standard library only): `python verify_outdegree_two.py`, `python verify_hardness_gadgets.py`, `python audit_hardness.py`, `python verify_variable_basis.py`, and `python verify_private_basis_spares.py`.
