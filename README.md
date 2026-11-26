# Robust Structural Controllability under Actuator Failures

Research on dedicated actuator placement under failures, including
approximation algorithms, exact algorithms, and fault-tolerant matroid
spanning. The state graph defines a transversal matroid; its dual translates
matching feasibility into a spanning condition. Source strongly connected
components impose additional accessibility quotas.

## Results and models

The [integrated arXiv manuscript](paper/robust-dedicated-actuator-placement.pdf) presents the 3/2 approximation, deterministic weighted source-quota algorithm, and one-sided randomized mixed-failure algorithm around the matching–accessibility decomposition.

- The [3/2 approximation](paper/robust-dedicated-actuator-placement.md)
  runs in cubic time on arbitrary directed state graphs for unit-cost
  dedicated actuators at distinct states under one actuator failure.
- [Factor-two approximation](two-approximation.md): for a nonempty directed
  state graph, unit-cost dedicated actuators at distinct states, and one
  actuator failure, a quota-first algorithm gives a factor-two guarantee and
  runs in cubic graph operations. The note also gives a family where the
  independent-repair strategy approaches factor two.
- [Bounded-core algorithm](bounded-core-theorem.md): the same one-failure
  placement model has an exact algorithm polynomial for each fixed matching
  deficiency `d`, obtained by enumerating robust cores of at most `2d` states.
- [Disjoint-quota theorem](quota-core-proof.md): unit-cost robust matroid
  spanning with disjoint quotas is fixed-parameter tractable in matroid rank
  and failure budget. Applied to dedicated actuators, the parameters are
  matching deficiency and the number of failures.
- [Weighted fault-tolerant spanning](paper/weighted-fault-tolerant-spanning.pdf)
  ([LaTeX source](paper/weighted-fault-tolerant-spanning.tex),
  [proof note](weighted-fault-tolerant-fpt.md)): weighted selected-set
  spanning in oracle matroids is fixed-parameter tractable in rank and
  failure budget, including the stated safe/vulnerable extension.

The [flat-saturation algorithm](weighted-flat-core.md) uses the quartic-time framework of Kundu, Nandi, and Saurabh and gives a cubic parameter exponent in cached independence queries, with subquartic regimes for restricted cardinality and packing. Bounds for special
classes appear in the [parallel-class theorem](parallel-class-fpt.md) and the
[partition-truncation theorem](truncated-partition-mixed-failures.md).

The [weighted purchase-quota theorem](weighted-disjoint-quotas.md) handles any
number of disjoint groups. The [robust-accessibility theorem](mixed-accessibility-quotas.md)
handles mixed failure statuses with an additional parameter for the source
groups. The [graph-resource family](graph-resource-approximation/README.md) has
an \(8/7\)-approximation. The [lower-bound excess theorem](lower-bound-excess-hardness.md)
gives NP-completeness at matching deficiency plus one.

The [linear source-quota algorithm](linear-source-quotas.md) ([PDF](paper/linear-source-quotas.pdf)) gives a
one-sided randomized exact algorithm for weighted actuator placement on
arbitrary directed state graphs, parameterized by matching deficiency
and failure budget.

The starting reference for robust structural controllability is Diogo Poças,
Sérgio Pequito, and Guilherme Ramos, [*On the NP-hardness of the robust
minimal structural controllability problem*](https://doi.org/10.1016/j.automatica.2026.113074),
*Automatica* 190 (2026), 113074.

## Sparse degree and hardness results

[Degree threshold and sparse hardness](hardness-degree-threshold.md): exact optimization for state outdegree at most two; APX-hardness with degrees in {2,3} and no mandatory actuator, including strongly connected bipartite graphs with undirected girth at least N^(1−delta) for every fixed 0<delta<1. The reduction preserves OPT=2|E|+tau and gives a 2238/2233−epsilon inapproximability bound using the absolute published Vertex Cover gap. [Reduction obstructions](hardness-obstructions.md) explain the failure of one-leaf compression, strengthen the forced-basis covering barrier to 121/108, and identify the published Activation Edge-Cover limitation of the smaller source-SCC baseline.

[Source-SCC covering reduction](source-scc-cover-hardness.md): pairing private basis elements in source two-cycles changes the cost relation to `OPT=2*(|U|+tau)` and gives **250/249−epsilon** hardness for general directed graphs. Source-loop completion removes all mandatory locations; exact odd padding gives undirected girth at least N^(1−delta). The targets remain loopless, bipartite and subcubic, with maximum indegree and outdegree three. All added states remain selectable. Unlike the preceding class, these graphs have several source SCCs; the two-copy construction preserves both solution directions without parity exceptions.

[Shared backup across a matching](shared-backup-hardness.md) improves the bounds to **913/908−ε** on general digraphs and **1594/1589−ε** on strongly connected digraphs. One resource can back up all jobs in a designated matching. Disjoint triangles use source quotas in the general construction and single matching jobs in the strongly connected one. Both reductions preserve the optional subcubic, bipartite and large-girth restrictions.

[Direct parity reduction](direct-lin2-hardness.md) gives **843/838−ε** hardness on general digraphs and **1474/1469−ε** on strongly connected digraphs with the same restrictions. Sharing a backup across literal-choice jobs leaves the forced resources below full matching rank. The excess above the explicit baseline equals the minimum number of violated Max-E3-Lin-2 equations.

[Stratified parity sampling](sampled-parity-hardness.md) gives **263/260−ε** on general digraphs and **503/500−ε** on strongly connected digraphs under **NP ≠ RP**. Three draws per variable and a bound on repeated hits of deleted variables give deletion cost above \(3n/20\) in the bad case. The target excess equals the minimum deletion count. An [anchored parity obstruction](parity-gadget-obstruction.md) proves that pure cotransversal matching needs at least seven fixed elements to encode ternary parity while independently prohibiting missing label pairs.

[Guard compression obstructions](guard-compression-obstructions.md) show when optional guards make all parity instances attain the matching lower bound. A local block with selectable auxiliary resources needs at least four jobs per independent ternary equation if failure of the shared backup pins all canonical labels. Run `python scripts/verify_guard_compression.py` for exact checks.

[General directed degree universality](general-sparse-equivalence.md): every feasible general directed instance, including multiple source SCCs and zero matching deficiency, has an exact realization in a weakly connected, loopless, bipartite graph with all underlying degrees in {2,3}, no mandatory location and girth at least N^(1−delta). Both solution maps are polynomial and cost nonincreasing. A connector made entirely of matching-matroid loops preserves accessibility without changing the optimum. More generally this realizes arbitrary coloop-free cotransversal robustness with disjoint quota-two groups that need not cover the ground set. The sparse class therefore has the same constant approximation threshold as the general problem; this is a structural theorem, not a new numerical lower bound.

[Partition-constrained assignment transfer](partition-assignment-transfer.md): arbitrary feasible unit-cost robust assignments with a partition of resource groups, each requiring two selections, have exact realizations in the same optional, bipartite subcubic control class with N-dependent girth and no added actuator baseline. One loop filler per resource makes source cycles bipartite without duplicating the instance. Reverse normalization retains selected original resources and hence all group quotas, allowing genuine choices beyond a fixed private basis.

[Clause reduction with rank-three matching blocks](clause-quota-hardness.md) proves NP-completeness at the accessibility lower bound of two actuators per source SCC, even when every component of the dual state matching matroid has rank at most three. Seven-element clause blocks and disjoint quotas encode satisfiability with variable literal choices. For bounded occurrence there is a quantitative gap transfer, with a smaller numerical bound than 250/249. A family of long variable cycles exhibits a single unit of excess cost even when the minimum violated fraction tends to 1/8.

[Clause gap transfer with bicircular copies](expander-clause-hardness.md) preserves a constant gap from OR or OXR formulas with unbounded variable occurrence. Minimum robust sets in small arity-two assignment blocks synchronize several ports. Directed expanders enforce consistency, and disjoint source quotas encode the clauses directly. The OXR reduction gives 2737/2736−ε hardness with dual matching components of rank at most twelve.

[Variable-basis path-cover equivalence](variable-basis-path-cover.md): an unforced construction with exact optimum `n+pc(H)`, polynomial solution maps, degree-three compression and the same girth padding. The classical (1,2)-TSP 8/7 algorithm limits this particular route to stronger hardness. [Private-basis hyperforest objective](private-basis-hyperforest.md) extends the construction to larger neighborhoods, proves `OPT=n+min_C(|C|-r_C)`, and explains why frequency-two and regular incidence sources lose their covering gaps.

[Sparse cotransversal transfer](cotransversal-hardness-transfer.md): every coloop-free positive-rank cotransversal matroid with a full-rank presentation of row size at most three has a strongly connected subcubic control realization with exactly the same robust optimum, polynomial solution maps and the same girth padding. A certified example is nontransversal for every basis; an explicit strongly connected source family refutes a head/tail-gate reduction from directed path cover. An acyclic gammoid-completion family has the sharp rank-preserving distortion 2r/(r+1), even though every added resource is individually avoidable in an optimum.

[Strongly connected degree universality](row-splitting-equivalence.md) removes the row-size restriction: every strongly connected directed instance has a loopless, bipartite, bidirected realization with degrees in {2,3}, no mandatory actuator, exactly the same OPT, cost-nonincreasing polynomial solution maps and the same girth padding. The underlying exchange applies to every coloop-free positive-rank cotransversal matroid. [Uniform assignment transfer](assignment-hardness-transfer.md) additionally handles arbitrary unit-cost uniform robust assignments without an actuator baseline; a dominated-resource exchange shows that adding a private basis preserves OPT and extends the hyperforest formula to sources without private resources.

[Basis-relaxation obstruction](basis-relaxation-gap.md): replacing a matroid by the assignment model of its fundamental exchanges at an arbitrary basis can underestimate robust OPT by a factor approaching 3/2, even for simple, cosimple fundamental transversal matroids with resources serving at most two jobs. Maximizing total column degrees in a maximal dual presentation does not fix the loss. This is a proved relaxation gap, not an NP-hardness claim. All 513 bases of the nontransversal prototype are also audited; every basis contained in an optimum placement has a strict relaxation gap in that example. A separate self-dual transversal family shows that every basis can admit spurious optimum-size placements even when all optimum values are exact.

Exact checks (Python 3.10+, standard library only): `python scripts/verify_outdegree_two.py`, `python scripts/verify_hardness_gadgets.py`, `python scripts/audit_hardness.py`, `python scripts/verify_variable_basis.py`, `python scripts/verify_private_basis_spares.py`, `python scripts/verify_gammoid_obstructions.py`, `python scripts/verify_assignment_transfer.py`, `python scripts/verify_row_splitting.py`, `python scripts/verify_basis_relaxation.py`, `python scripts/verify_scc_cover.py`, `python scripts/verify_partition_assignment.py`, `python scripts/verify_quota_cotransversal.py`, `python scripts/verify_clause_quota.py`, `python scripts/verify_amplified_clause_quota.py`, `python scripts/verify_shared_backup_cover.py`, `python scripts/verify_direct_lin2.py`, and `python scripts/verify_parity_guards.py`.

## Verification

Python 3.10 or newer and the standard library are sufficient. From the
repository root:

```sh
python scripts/verify_core_theorem.py
python scripts/check_two_approx.py
python scripts/audit_weighted_core.py
python scripts/audit_weighted_independent.py
```

The checks compare algorithms with exhaustive or independent finite-instance
oracles. Each proof note specifies its test domain and corresponding solver.
The [paper guide](paper/README.md) gives manuscript build instructions and
links to the integrated control manuscript.
