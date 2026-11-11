# A 3/2 approximation via maximum co-repair matching

## Model and baseline

The state graph is nonempty and directed, with self-loops allowed. Costs are one, actuators are dedicated, each selected state receives at most one actuator, and any one actuator may fail. Failure removes the selected input only, not a state vertex or any state edge. Write E for the state set and n=|E|. A feasible instance has no singleton source SCC. Let T be the transversal matroid of right-state subsets matchable by state edges, M=T*, and d=rank(M). A placement is robust precisely when it remains spanning in M after any actuator deletion and meets every source SCC at least twice.

Use the exact quota-spanning relaxation in [the factor-two note](two-approximation.md): minimize |P| subject to spanning M and placing at least two actuators in every source SCC. Its max-flow construction gives L=|P|≤OPT and a matching covering E minus P. Complete repair supports, explored in [PR #2](https://github.com/dongxuelian2/robust-structural-controllability/pull/2), are used below. This implementation depends only on merged baseline files, not on PR #2 being merged. Publication novelty remains unverified.

## Co-repair graph

Let M be a matroid of rank d, P a spanning set, and D its restricted coloops. For b in D put C_b=E\cl(P-b). For e outside P define A_e={b in D:e in C_b}. Define G on D with edge bc whenever some A_e contains both b and c. Let nu be the maximum cardinality of a matching of G.

## Rank identity for arbitrary R subset D
Choose a basis B subset P. Every b in D is in B. Since P-b and B-b both have rank d-1 and B-b is contained in P-b, their closures agree. Thus C_b is the fundamental cocircuit of b with respect to B.

We claim intersection over b in R of cl(B-b) equals cl(B-R). The inclusion from right to left is immediate. For an element e outside B, membership in cl(B-b) is equivalent to b not belonging to the unique fundamental circuit of e in B+e. Satisfying this for every b in R says that this circuit is contained in (B-R)+e, equivalent to membership in cl(B-R). For e in B, both sides contain exactly B-R. (Loops are contained in every closure.) Therefore

rank(E minus union_{b in R} C_b)=d-|R|.

This proves the delicate rank equality without assuming representability or modularity of arbitrary matroid flats.

## Component lower bound
Fix any U subset D and let R_1,...,R_s be the nonempty connected components of G-U. Put Q_i=union_{b in R_i}C_b. These Q_i are disjoint: a common external element would cover b and c from two different components and create an edge; a point of P belongs to only its own C_b.

Let W be any single-failure robust spanning set. Since E\Q_i has rank d-|R_i|, W must contain at least |R_i| elements of Q_i. If it contained exactly |R_i|, deletion of any one of them would leave rank at most (d-|R_i|)+(|R_i|-1)<d. Hence

|W intersect Q_i| >= |R_i|+1.

Summing over disjoint Q_i gives |W|>=|D|-|U|+s, which is at least |D|-|U|+odd(G-U).

The Tutte–Berge formula gives

2nu=|D|-max_U(odd(G-U)-|U|).

Taking the maximum of the preceding lower bound yields |W|>=2(|D|-nu).

## Repair and guarantee
Choose a maximum matching in G. For each matching edge choose one common external repair; for each unmatched b choose any external repair. Their union R repairs all of D and uses at most |D|-nu elements (duplicated choices only reduce size). Thus |R|<=|W|/2 for every robust W, including an optimum placement.

In the controllability setting, choose P as the exact minimum quota-spanning relaxation, L=|P|<=OPT. Adding R preserves the source-SCC quotas, and complete coverage makes the spanning property robust. Hence |P union R|<=L+OPT/2<=3OPT/2.

The same argument also proves that *exact* minimum repair of any such fixed P costs at most OPT/2. Exact set cover need not be solved: maximum matching supplies this upper bound constructively in polynomial time. All statements use ordinary unit cardinalities.

## Tight example
The three-vertex digraph with outgoing masks [6,1,1] (arcs 0->1,2 and 1,2->0) gives an exact 3/2 example for the current quota-flow tie breaking: P={0,2}, D={2}, repair={1}, output size3, optimum {1,2} of size2. Thus arbitrary optimum-relaxation choice plus this repair approach cannot promise a factor below3/2. Reselecting P can solve this example, so it is not an inapproximability statement.

## Polynomial implementation and reproducible checks

`three_halves_solver.py` retains the exact quota-flow matching from `two_approx_solver.py`. For each critical b, alternating search from unmatched right b visits exactly the external states e that can replace b. A reachable free left endpoint would contradict criticality. A path flip proves sufficiency; the symmetric difference with a replacement matching proves necessity. Creating a clique on each repair support constructs the co-repair graph and records a common repair witness for each edge. An Edmonds augmenting-forest algorithm with blossom contraction computes a **maximum**, not merely maximal, matching. Matched pairs and unmatched vertices then select repair witnesses, deduplicated before installation.

The quota flow, at most n alternating searches, clique construction and general-graph matching take O(n³) graph operations in total, using O(n²) storage. Python integer bitmasks have their usual bit-complexity costs. A maximal matching or a greedy selection of pairs does not justify the Tutte–Berge lower bound used here.

Use Python 3.10+, without third-party packages:
- `python audit_blossom.py`: compares actual blossom matching with independent subset-DP maximum matching on every simple graph of orders 0 through 6 and 1,200 seeded random graphs of orders 7 through 18, 35,068 total.
- `python audit_three_halves.py`: all 66,066 looped state graphs through order 4 and 3,500 seeded random graphs of orders 5 through 11, 69,566 total. It independently enumerates matchable subsets, robust optima and repair supports, checks the maximum co-repair matching against subset DP, and verifies returned matching/repair certificates and the 3/2 bound.

- `python audit_three_halves_independent.py`: 66,066 exhaustive plus 800 random instances; independently checks every optimum P, every critical subset rank identity, every Tutte–Berge vertex deletion, and the actual blossom solver against subset DP. No unpublished PR #2 helper is required.

These finite checks supplement the proof; they do not establish publication novelty.

## Related work and limits

Adjiashvili, Bindewald and Michaels, *Robust Assignments with Vulnerable Nodes*, [Theorem 7](https://arxiv.org/pdf/1703.06074), proves 7/4 for uniform cardinality node-robust assignment. Its greedy productive/singleton analysis informed the preceding joint-repair work; the present proof uses maximum co-repair matching and Tutte–Berge instead.

The same authors' *Robust Assignments via Ear Decompositions and Randomized Rounding*, [ICALP 2016, Theorem 5](https://doi.org/10.4230/LIPIcs.ICALP.2016.71), already proves 3/2 for uniform cardinality **edge**-robust assignment. That model buys assignment edges and tolerates an edge failure, using ear decompositions. It is not the dedicated-actuator/node-failure problem here, and that theorem is not applied as a black box.

Other relevant sources are [Poças, Pequito and Ramos, Automatica 190 (2026), 113074](https://doi.org/10.1016/j.automatica.2026.113074) for robust structural controllability, and [Bentert, Fomin, Golovach and Morelle, Fault-Tolerant Matroid Bases, ESA 2025](https://doi.org/10.4230/LIPIcs.ESA.2025.83) for the abstract robust-spanning setting. Maximum matching, fundamental matroid circuits/cocircuits and the Tutte–Berge formula are established tools. No first-result or best-known claim is made; novelty of this particular combination remains unverified.

The result makes no weighted-cost, multiple-failure, unrestricted-input or repeated-actuator claim. The 3/2 equality example limits this relaxation-and-repair choice, not every algorithm for the underlying problem. Exact minimum fixed-P repair is not computed or required.
