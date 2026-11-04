# Subcubic bidirected graphs preserve the strongly connected problem

**Theorem.** Every strongly connected directed instance on at least two
states has a polynomial-size, loopless, bipartite, bidirected, strongly
connected instance of maximum degree three with exactly the same robust
dedicated-actuator optimum. Both solution maps are polynomial and do not
increase cost. For every fixed 0<delta<1, the target can additionally have
undirected girth at least N^(1-delta). All actuators have unit cost, remain
distinct, and are vulnerable to failure.

More generally the same exact reduction applies to robust spanning in
**every coloop-free positive-rank cotransversal matroid**, given any
transversal presentation of its dual. There is no row-size restriction and
no added actuator baseline. With q presentation incidences, the unpadded
control target has fewer than 10q states.

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
cost nonincreasing. Odd-subdivision girth padding preserves these maps.
If T is simple, each split T' remains simple: original independent pairs
stay independent under a coextension, and x has two row neighbors, so it
is independent of every nonloop original element. Hence the construction
also preserves the absence of two-element cocircuits in M, and therefore
the absence of mandatory actuators in the final control graph.

For a strongly connected control source on at least two states, take T to
be its state transversal matroid. If the matching deficiency d>0, a robust
spanning set of T* has at least d+1>=2 elements; the sole source-SCC quota
adds no further condition. Also T has no loops. The theorem follows from
the preceding construction, with q the number of source arcs, including
self-loops. If d=0, every pair of states is robust, so OPT=2; use a
bidirected two-state edge as the target and map its two actuators to any
two source states. Its undirected girth is infinite.

This statement concerns strongly connected graphs. Multiple source-SCC
accessibility quotas in general directed graphs require a separate
reduction and are not absorbed into the matroid argument.

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

`cotransversal_reduction` and `cotransversal_from_placement` in
`row_splitting.py` implement the two maps. The reverse map also accepts
instances padded by `subdivide_ports`.
