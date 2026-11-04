# Sparse cotransversal matroids as control instances

The private-basis formula in `private-basis-hyperforest.md` concerns a
fundamental transversal matroid. This note gives a reduction for a larger
class, a certificate that the larger class is actually needed, and an
obstruction to a proposed directed-path reduction. It establishes no new
numerical inapproximability constant.

The incidence construction below assumes row size at most three.
`row-splitting-equivalence.md` now removes that assumption by an exact
row-splitting exchange and handles arbitrary, initially rank-deficient
presentations. In particular every strongly connected control instance
has an exact subcubic bidirected realization.

## A general reduction with no actuator baseline

Let T be a transversal matroid on n elements, given by r full-rank rows,
each containing at most three elements. Suppose T has no loops and 0<r<n.
Then M=T* is coloop-free with positive rank d=n-r. Write rho(M) for the
minimum size of a set that spans M after any one selected-element deletion.

**Theorem.** In polynomial time construct a loopless, bidirected, bipartite,
strongly connected control graph G of maximum degree three such that

    OPT(G) = rho(M).

Both solution maps are polynomial and do not increase cost. The initial
construction has at most 11r-n states. For every fixed 0<delta<1 it can also
have undirected girth at least N^(1-delta), still with exactly the same
optimum and solution maps. If M has no two-element cocircuit, no actuator
in G is mandatory. A planar presentation can be realized planarly when its
embedding orders are supplied.

**Proof.** First use the bidirected incidence graph of the presentation,
with n element targets and r row centers. Its state transversal matroid is

    T on targets  direct sum  a free matroid on centers.

Indeed left centers alone match target columns according to T. A full row
matching in the presentation, transposed, matches every right center using
distinct left targets. The two row families are disjoint, so every target
independent set can be extended by all centers. Thus the dual matching
matroid is M together with r loops; the deficiency is d.

Split each degree-k element target into a path of k copies separated by
k-1 new degree-two centers, distributing one original incidence to each
copy. Each target copy has degree at most three; original centers already
have degree at most three. The target-splitting lemma in
`variable-basis-path-cover.md` shows that each added target is a parallel
copy in the dual matching matroid, and each new center is a loop. This is
a presentation operation for arbitrary T, independent of a private basis.

Each presentation component has a center of degree at most two after
splitting, unless it was just one degree-three center with three leaves.
In that exceptional case, subdivide one incidence by a length-three path
to create a degree-two center. At a chosen center port, attach a new target
q and a new center p by port--q--p. In the target presentation p has only
neighbor q, so q is a coloop of the transversal matroid; both new states
are loops of the dual matching matroid. Join the p centers of successive
components by bridges. Their degrees stay at most three. These additional
center--center edges do not affect the freely matchable center columns.
The graph is connected and bidirected, hence strongly connected. Component
bipartitions can be flipped successively along the bridge path.

Let E be the original incidence count and c the number of components.
After target splitting there are r+2E-n states and 3E-2n edges. Ports and
bridges add at most 4c states and 5c-1 edges. Since E<=3r and c<=r,

    N <= 11r-n,       |E(G)| <= 14r-2n-1 < 12n.

The resulting matching matroid is a parallel extension of M plus loops.
The coloop-completion lemma in `hardness-obstructions.md` maps any robust
selected copies back to a robust original set of no larger size. Original
robust sets remain robust using their retained representatives. Every
robust set has at least d+1>=2 elements; accessibility in the connected
control graph therefore adds no cost. This proves the optimum identity
and the reverse map, implemented by `protect_coloops`. The constructor is
`cotransversal_graph` in `hardness_gadgets.py`.

Parallel extension expands cocircuits by replacing members with their
parallel classes; it cannot create a two-element cocircuit when none
existed. Added loops are never mandatory. This proves the assertion about
mandatory actuators. Target splitting follows the cyclic incidence orders
of a planar embedding, and component bridges preserve planarity.

Finally subdivide every center--target edge to length L=3^k, leaving only
the center--center connector bridges unsubdivided. All cycles stay within
presentation components. Consequently N<=12nL and girth>=4L. For a fixed
integer a>=max(1,(1-delta)/delta), take the least such L>=(12n)^a. Then N is
polynomial in n and 4L>=N^(1-delta). The same parallel normalization applies
after any number of rounds. This completes the proof.

This transfers a hardness gap on the stated cotransversal class with no
dilution: it is an L-reduction with constants (1,1). It does not assert
that every strict gammoid has a presentation with row size at most three.
The cotransversal/strict-gammoid duality itself is classical; see
[Fink--Olarte, *Presentations of transversal valuated matroids*, Section 7](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/jlms.12505).
No literature-priority claim is made for the incidence realization.

## A connected example outside every private-basis representation

On elements 0,...,12 take the following eight rows of T:

    {0,3,10}, {2,9,12}, {2,5,11}, {3,5,10},
    {4,5,9}, {5,6,8}, {6,7,8}, {1,4,12}.

They have full row rank eight. Their dual M has rank five, is connected,
has no parallel elements and no two-element cocircuits, and has robust
optimum seven. A placement of size seven is

    {0,1,2,3,4,6,7}.

The finite certificate in `verify_gammoid_obstructions.py` enumerates all
8,192 supports, checks all 1,716 six-input exclusions independently using
actual state-plus-input matching, and checks every failure of this witness.
The 21-state incidence graph has one degree-four target. Splitting target
5 between rows {2,3} and {4,5} adds just two states and yields a subcubic
23-state graph with the same rank and optimum. The script verifies the
complete 16,384-entry parallel-extension rank identity.

More significantly, M is **not transversal**, regardless of basis choice.
Consider these cyclic flats:

    A={5,6,7,8},   B={1,4,12},   C={2,9,12}.

Their ranks, computed directly by the dual matching formula, are:

| Set | A | B | C | A union B | A union C | B union C | A union B union C | A intersect B intersect C |
| --- | - | - | - | --------- | --------- | --------- | ----------------- | ------------------------- |
| Rank | 2 | 2 | 2 | 4 | 4 | 3 | 4 | 0 |

For a transversal matroid, Mason's cyclic-flat inequality requires

    r(A intersect B intersect C)
      <= r(A)+r(B)+r(C)-r(A union B)-r(A union C)-r(B union C)
         +r(A union B union C).

Here it would require 0<=-1. The necessary inequality, including its
cyclic-flat hypothesis, is [Bonin--Kung--de Mier, *Characterizations of
Transversal and Fundamental Transversal Matroids*, Theorem 3.2](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v18i1p106/pdf/).
Thus the example lies outside the entire private-basis assignment family,
not just outside the representation for a particular basis. Its sparse
control realization remains nontransversal: deleting the new parallel
copies and loops recovers M, and transversality is closed under deletion.
This is a test instance for genuinely cotransversal reductions, not an
asymptotic approximation lower bound.

## Why head and tail gates do not preserve directed path cover

Here is a different candidate, and a linear counterexample to its proposed
cost map. Given a directed source H, make sink states b_v for its vertices
and two more sinks h_0,h_1. Use gate vertices

    f_v -> {b_v,h_0},       g_v -> {b_v,h_1},

and for every source arc u->v an actuator candidate a_uv with arcs

    a_uv -> {f_u,g_v}.

Take the strict gammoid linking all selected endpoints to the sinks. Gates
are selectable like every other element. Its rank is |V(H)|+2. A hoped-for
map OPT=|V(H)|+pc_directed(H)+O(1) would transfer a directed TSP gap with
only constant overhead. It is false even for strongly connected sources.
Failure removes a starting resource; its graph vertex remains available
for other paths. The dual transversal presentation has one row
{w} union N^+(w) per non-sink w. Here every row has size three and the
identity entries give full row rank, so the theorem above turns these
networks into control instances with exactly the same optimum.

For k>=4 take a directed backbone v_0->...->v_(k-1), two leaf children
l_i,r_i of every v_i, and a return vertex z. Add leaf->z for every leaf
and z->v_0. This source is strongly connected, has 3k+1 vertices, and

    pc_directed(H)=2k-1.

Indeed all 2k leaves have only z as an out-neighbor, so at most one can be
internal in a path partition. The bound is attained by the path
l_0,z,v_0,...,v_(k-1),l_(k-1), with every other leaf a singleton.

Yet the target's robust optimum is exactly 4k. For the lower bound, the
only starting resources able to reach a leaf sink b_l are

    {b_l,f_l,g_l,a_(parent,l),a_(l,z)}.

These five-element sets are disjoint for the 2k leaves. Each needs at least
two selected resources, since its only selected resource could fail.

For the upper bound, for each branch select a_(v_i,l_i), a_(v_i,r_i),
f_(l_i), f_(r_i). Replace f_(l_0) and f_(l_1) by a_(l_0,z) and a_(l_1,z).
There are 4k selected elements. After any failure, use one of the replacement
arcs in a different branch to cover b_z. The failed branch and the branch
serving z each retain three resources, enough for their three local sinks.
Every intact four-resource branch can cover its three local sinks and
either h_0 or h_1. At least two branches remain intact, so assign the hubs
to two of them. These paths use disjoint local gates; only their distinct
hub endpoints are shared across branches. This gives a full linking after
every failure. Replacement arcs not assigned to z follow their old f_l
route. Therefore

    OPT=4k,       |V(H)|+pc_directed(H)=5k.

The discrepancy grows with k and cannot be absorbed into an additive
constant. Actuators on internal gates undo the intended orientation.

## Verification and next requirement

`python verify_gammoid_obstructions.py` is standard-library only. In
addition to the finite nontransversality certificate, it checks the comb
linkings by independent integer vertex-capacitated max flow for k=4,...,12.
It exhausts all 1,595 loopless full-row-rank presentations on at most four
elements with row size at most three, including disconnected presentations
and zero-dual-rank components. Four further padded fixtures give 1,599
instances, 80,365 full-state rank queries and 15,150 reverse-map checks.

A stronger hardness proof still needs an asymptotic gap within this larger
cotransversal class. The realization theorem removes the graph-degree and
actuator-baseline obstacles; the remaining source construction must retain
its gap when every internal routing vertex can itself be selected.
