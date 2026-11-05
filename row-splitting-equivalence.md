# Subcubic bidirected graphs preserve the strongly connected problem

**Theorem.** Every strongly connected directed instance on at least two
states has a polynomial-size, loopless, bipartite, bidirected, strongly
connected instance with every degree in {2,3} and exactly the same robust
dedicated-actuator optimum. Both solution maps are polynomial and do not
increase cost. For every fixed 0<delta<1, the target can additionally have
undirected girth at least N^(1-delta). All actuators have unit cost, remain
distinct, and are vulnerable to failure. No actuator is mandatory.

More generally the same exact reduction applies to robust spanning in
**every coloop-free positive-rank cotransversal matroid**, given any
transversal presentation of its dual. There is no row-size restriction and
no added actuator baseline. With q presentation incidences, the unpadded
base control target has fewer than 10q states. Removing its leaves as
below gives the degree-{2,3} version with fewer than 50q states.

This removes the source restriction in `cotransversal-hardness-transfer.md`.
It is a structural equivalence, not a new numerical approximation gap.

## Compressing a presentation to full row rank

Let T be presented by any number p of rows on n target elements. Find a
maximum row-to-target matching, of cardinality r. Starting from all
unmatched rows, follow alternating paths: row-to-target edges and matched
target-to-row edges. Write R_0,V_0 for the reached rows and targets. All
targets in V_0 are matched, and their matched rows belong to R_0. Also
N(R_0)=V_0. The rows in R_0 therefore match every target in V_0
simultaneously, and cannot match any target outside V_0.

Replace R_0 by one singleton row {v} per v in V_0. This preserves every
matchable target subset: on V_0 these rows supply a free matroid, and the
other rows are untouched. The number of rows is now

    p-|R_0|+|V_0| = p-(p-r) = r.

The new presentation has full row rank and no more incidences than before.
This elementary alternating-matching normalization is implemented by
`full_rank_presentation` in `row_splitting.py`.

## Splitting a large row preserves the robust optimum

Let one presentation row R be the union of two nonempty sets A,B. Replace
it by rows A union {x} and B union {x}, where x is a new target element.
Let T' be the resulting transversal matroid and M'=T'^*. Then

    T'/x = T,             M' restricted to the old ground set = M=T*.

Indeed a matching required to include x uses one new row for x and the
other for whichever element formerly used R. Conversely contracting x
leaves at most one original target assigned to the two split rows. Thus
r(T')=r(T)+1 and r(M')=r(M). Full row rank and the absence of target loops
are preserved.

**Row-splitting lemma.** Every robust spanning set of M' can be converted
to a robust spanning set of M of no larger size, in polynomial time.

**Proof.** For an actuator set S, put I=V\S. Matroid duality says S is
robust spanning iff I is an independent flat of T: I is independent, and
I+e is independent for each e outside I.

Consider a robust set S' of M'. If x is absent, the restriction identity
already gives an original robust set. Otherwise put S=S'\{x} and
I=V\S on the original ground set. Failure of x says S spans M, so I is
independent in T. For every e in S such that I+e is dependent in T, let
C_e be its unique circuit. Robustness of S' says I+e is independent in T',
so this circuit's Hall obstruction must use the changed row R. Minimal
dependence gives

    |N_T(C_e)|=|C_e|-1,
    N_T(C_e\{e})=N_T(C_e),
    R belongs to N_T(C_e\{e}).

Thus K_e=C_e\{e} is a tight Hall subset of the independent set I, and
every N_T(K_e) contains R. Tight Hall subsets of I are closed under union
and intersection, with

    N_T(K_1 intersect K_2)=N_T(K_1) intersect N_T(K_2).

This follows from neighborhood submodularity and Hall on every subset of
I. Consequently the intersection of all K_e is nonempty: its neighborhood
still contains R. Choose w in that intersection. Removing w destroys every
formerly dependent I+e circuit. Every other I+e was already independent,
and (I\{w})+w=I is independent. Therefore I\{w} is an independent flat
of T, giving the robust original set S union {w} of size |S'|. If there
are no such circuits, S itself is already robust. To implement the
exchange, try deleting x, then try each original unselected w, testing
robustness by matching. No circuit or intersection oracle is needed.

Conversely an original robust set stays robust under the restriction
identity. This proves equality of the two optima and both solution maps.

The lemma also explains why a new branch vertex cannot silently become a
cheaper actuator. Its selected unit can always be exchanged back to one
original element; no assumption about mandatory actuators is used.

The exchange has a stronger local form. In any presentation of a loopless
transversal T, a column x with at least two row neighbors can be deleted
from M=T* without changing its robust optimum. To prove this, start with
robust S containing x and put I=V\S, J=I+x. The set J is independent.
For each e in S\{x} with J+e dependent, its unique circuit C_e contains
x, since I+e is independent. The sets C_e\{e} are tight Hall subsets of
J containing x. Their intersection K is tight and contains x. Because
N(K) contains both row neighbors of x, |K|>=2. Any w in K\{x} removes
every obstructing circuit, so S\{x} union {w} is robust on the old ground
without x. If no circuit obstructs, simply delete x. The converse follows
by restriction. `avoid_presentation_element` implements this exchange.
The resulting deletion need not itself be cotransversal.

This is an individual deletion statement. It cannot be iterated to avoid
an arbitrary prescribed set: the completion family in
`cotransversal-hardness-transfer.md` has optimum r+1 although forbidding
all added vertices raises it to 2r, and every added column individually
has at least two row neighbors.

## Removing leaves without changing the optimum

In any of the incidence realizations, attach a new four-cycle by one
bridge at every old leaf v. The cycle has two center states and two
target states. If v is a target, join it to a cycle center; if v is a
center, join it to a cycle target. Every old leaf now has degree two, the
cycle attachment point has degree three, and its other vertices have
degree two. No new leaf is introduced. Bipartiteness, connectivity and
planarity are preserved. This is `complete_leaves` in `hardness_gadgets.py`.

Right centers remain a free summand: new left targets perfectly match
the new right centers, separately from target matching. If v was a center,
the two new target columns form a free summand as well. The two new row
centers see only the new targets, so they match both without consuming
any original row. All four added states are therefore dual loops.

If v was a target, the two new targets occur only in the two new rows.
Contracting them in the target transversal presentation recovers the old
presentation, so the original dual matroid is a restriction of the new one.
Deleting v and either one of the new targets leaves only one new target
available to the two new rows. The dual rank formula then gives the same
rank for this pair as for v alone. If v is a dual nonloop, a full original
row matching avoiding v shows that either new target is also a nonloop.
Thus both are parallel copies of v. If v was a dual loop, both are loops.
The added centers remain dual loops.

Consequently the entire operation adds only parallel copies and loops.
Coloop completion gives a polynomial reverse map of no larger cost.
It adds four states and five edges per old leaf, hence at most multiplies
the state count by five. With minimum degree two the state transversal
matroid is simple: every pair of target columns has a two-row matching.
Equivalently, on the strongly connected target, actuating every state
except any prescribed v is robust after any other failure. This proves
that no actuator is mandatory, including when the source had forced ones.

This removes mandatory locations, while preserving their cost through
the parallel-copy normalization. It does not bypass the forced-basis
approximation barrier for a source that already has that architecture.

## General cotransversal and control consequences

For a row R of size k>3 split off two elements into A, let B=R\A, and
apply the operation. One new row has size three and the other size k-1.
Repeat until all rows have size at most three. For a full-rank original
presentation with r rows, n elements and q incidences this uses

    a=sum_R max(0,|R|-3) <= q-r

splits. The resulting presentation has r+a rows and n+a elements, with
unchanged dual rank and unchanged robust optimum. Apply the subcubic
incidence realization from `cotransversal-hardness-transfer.md`. Its state
count satisfies

    N <= 11(r+a)-(n+a) <= 10q+r-n < 10q.

Its parallel copies and loops are removed first by coloop completion;
then reverse the row splits using the lemma. All maps are polynomial and
cost nonincreasing. Apply the leaf completion to obtain N<50q and at most
62q edges. Odd-subdivision girth padding preserves these maps and the
degree-{2,3} condition. With length multiplier L=3^k, N<=62qL and
girth>=4L. Fix an integer a>=max(1,(1-delta)/delta) and choose the least
such L>=(62q)^a; the same inequality as in the sparse incidence theorem
gives girth>=N^(1-delta), with polynomial graph size.
If T is simple, each split T' remains simple: original independent pairs
stay independent under a coextension, and x has two row neighbors, so it
is independent of every nonloop original element. Hence the construction
also preserves the absence of two-element cocircuits in M before the leaf
completion. The completion removes mandatory actuators without this
additional hypothesis.

For a strongly connected control source on at least two states, take T to
be its state transversal matroid. If the matching deficiency d>0, a robust
spanning set of T* has at least d+1>=2 elements; the sole source-SCC quota
adds no further condition. Also T has no loops. The theorem follows from
the preceding construction, with q the number of source arcs, including
self-loops. If d=0, every pair of states is robust, so OPT=2; use a
bidirected even cycle on at least four states as the target and map any
two target actuators to any two source states. Its girth equals its order,
and no target state is mandatory.

This statement concerns strongly connected graphs. The separate construction
in `general-sparse-equivalence.md` now retains multiple source-SCC quotas
as well, using new source representatives and a connector of dual loops.
Its target is weakly connected and need not be bidirected. It gives exact
degree-{2,3} universality for the general directed problem.

## Verification

`python verify_row_splitting.py` checks presentation compression on all
3,241 multisets of at most six rows on up to three targets, including zero
rows and rank-deficient presentations. It compares the entire transversal
rank table before and after compression. It then checks 298 row-splitting
sources, including all 132 feasible row multisets on four elements that
have a size-four row, eight contractions of the certified nontransversal
rank-five prototype, and larger fixed-seed examples. In total it checks
53,232 original-restriction ranks and 22,061 reverse maps. It separately
checks simplicity preservation, the actual input-failure model on sparse
and padded targets, and five strongly connected directed source graphs
whose original presentations have redundant rows.
The stronger individual-deletion exchange is checked on all 448 loopless
full-row-rank row multisets on up to four targets, with 9,046 exchanges.
Leaf completion is checked on 14 sparse or padded realizations, with
1,120 full-state rank queries and 1,032 reverse maps. The separate cover
gadget verifier checks eight leafless or padded realizations and 404
robust cover-extraction maps.

`cotransversal_reduction` and `cotransversal_from_placement` in
`row_splitting.py` implement the two maps. The reverse map also accepts
instances padded by `subdivide_ports`.
Compose `complete_leaves(cotransversal_reduction(...))` for the minimum
degree-two version; the same reverse map accepts this completion.
