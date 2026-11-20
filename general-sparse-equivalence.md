# General directed instances have exact sparse realizations

**Theorem.** Every feasible general directed instance of unit-cost,
dedicated, single-selected-actuator-failure controllability has a
polynomial-size realization that is weakly connected, loopless and
bipartite, with every underlying undirected degree in {2,3}, both directed
degree bounds three, and no mandatory actuator location. The two instances
have exactly the same optimum and polynomial cost-nonincreasing solution
maps. For every fixed 0<delta<1 the target can also have undirected girth
at least N^(1-delta). All added states remain selectable.
The matching deficiency and the number of source SCCs are preserved.

Thus the optimal constant approximation threshold of this sparse class
equals that of general directed graphs. This strengthens the previous
strongly connected universality theorem: multiple source-SCC quotas and
unconstrained nonsource states are retained. The target need not be
strongly connected or bidirected. This is an exact structural reduction,
not a new numerical hardness gap; 3/2-epsilon hardness remains open.

We prove the following more general statement. Let M=T* be any coloop-free
cotransversal matroid on E, |E|=n, given by a transversal presentation of
T. Impose quota two on disjoint groups P_1,...,P_t, where t>=1 and each
group has at least two elements. The groups need not cover E. Robust
spanning with these quotas has an exact realization in the class above.
Rank zero is allowed. If q is the number of presentation incidences, the
unpadded target has at most 13q states.

For a control source take its state transversal matroid T and its source
SCCs as the groups. A feasible instance has every source SCC of size at
least two. Every state then has an incoming arc, so T is loopless. Matroid
duality identifies the matching condition with robust spanning in T*;
accessibility under one failure is precisely quota two in each source
SCC. This also gives a polynomial feasibility check. The abstract theorem
therefore proves the control theorem, including matching deficiency zero.

## Binary row trees with available root ports

Compress the presentation to full row rank k without increasing its
incidence count, as proved in `row-splitting-equivalence.md`. Choose a
matching of these k rows to distinct original elements p_i. For a row
R_i of size s>=3 split it first into

    {p_i,x},             (R_i\{p_i}) union {x}.

Subsequent splits peel off two elements, including the current matched
column, until the remaining rows have size three. This uses s-2 splits.
Rows of size one or two are unchanged. The root row has size at most two,
each other row has size three, and a new column joins adjacent row-tree
nodes. Each split preserves the old matching and matches the appended
row to the new column. The root still matches p_i.

Writing a=sum_i max(0,|R_i|-2), the resulting T' has k+a rows and n+a
columns. It is a coextension of T, so M'=T'^* restricts to M and has the
same rank d=n-k. The row-splitting exchange supplies a cost-nonincreasing
reverse map. In particular it never removes a selected original element:
it removes only the last added column, and either stops or adds one old
element. For completeness, the exchange uses the independent-flat dual
description of robustness. The tight Hall subsets obstructing deletion of
the split element all meet the changed row; their intersection contains
an old element whose exchange removes every obstruction. The full proof
and a matching-based implementation are in the preceding note.

## An oriented incidence realization

Create a physical target v for each column of T', and a center c_i for
each row. Put c_i->v for each presentation incidence. For the retained
full row matching, also put p_i->c_i. Unmatched physical targets have no
outgoing arc. Color physical targets zero and centers one.

The center columns form a free transversal summand: the distinct matched
physical rows match them all. The rows at centers present exactly T' on
physical columns. Consequently the dual state matching matroid is M'
plus loops on the centers. Each root center has underlying degree at
most two, since its incoming matched arc reverses an existing incidence.
Each other center has degree three. Its matched target has just that
center as an outgoing neighbor, giving a private bidirected gate edge.
Every new row-tree center is reachable from its root through these gates.

## Source groups and a connector made entirely of matroid loops

For every grouped original element e make a source representative a_e
and a filler f_e. For an ordered group (e_1,...,e_s), use the cycle

    a_(e_1) -> f_(e_1) -> ... -> a_(e_s) -> f_(e_s) -> a_(e_1).

From the predecessor filler of a_e add an arc to physical target e.
Color a_e zero and f_e one. Ungrouped original elements retain their
physical targets as their original representatives.

After ignoring the free color-one columns, the added row is {a_e,e}.
This is a series extension of T', hence a parallel extension of M'.
Each a_e is a parallel copy of e (or a loop if e is a loop). The fillers
are matroid loops. These cycles impose exactly the desired group quotas.

Some physical targets are ungrouped, so accessibility requires an extra
connector. From one a_e in the first source cycle, reach all k root
centers using a binary tree. Each internal pair is x of color one and
y of color zero, with

    parent -> x,       x -> y,       y -> x,

and at most two child arcs leaving y. Leaves are the root centers.
There are 2(k-1) new states; for k=1 use one direct arc instead.

**Every connector state is a matching-matroid loop.** The y columns are
free, using their singleton x rows. All color-one columns are still
simultaneously matchable: physical matched rows serve the old centers,
source a rows serve the fillers, and connector y rows serve their x
partners. Additional arcs to centers affect only this already free
summand. Thus the connector adds only free transversal columns, or dual
loops, and changes no functional ranks. It must not be confused with a
gammoid routing tree that makes its internal vertices useful actuators.

The connector uses one spare underlying port at a_e and at each root
center, so degrees remain at most three. It reaches every physical target
and every center. The source SCCs are exactly the displayed cycles;
there are no other source SCCs. The underlying graph is connected even
if the presentation was disconnected. It is not necessary to force any
actuator on the connector.

## Completing the degrees

Only original physical targets can exceed degree three. A matched target
has one outgoing arc, reversing a retained incoming matching arc; protect
that incoming arc. While more than three incoming arcs remain, move two
other incoming arcs to a new x and add

    x -> y,       y -> x,       y -> v.

The new x is parallel to v in the matching matroid and y is a loop, by
the column-splitting lemma. The indegree of v decreases by one, and its
outgoing neighbor remains one of its incoming neighbors. All underlying
degrees are at most three. The operation preserves the bipartition and
the source SCCs.

The only possible leaves are ungrouped original physical targets. At
each such v attach a bidirected four-cycle by a bridge to a cycle center.
The two new targets are parallel copies of v and the new centers are
loops. More generally any even cycle with k targets and k centers works:
contracting its targets in the target transversal recovers T', while
deleting v and any one new target loses one row rank. If v is a dual
nonloop, match the bridge center to v and the other k-1 cycle rows to
the surviving cycle targets, with the old rows matched avoiding v.
Deleting any row and target of a cycle leaves this near-perfect matching.
Thus every new target is parallel to v; for a dual loop all are loops.
No new leaf is introduced.

At this point every incoming column neighborhood is nonempty. A column
with just one incoming neighbor is a source-cycle column or a private
gate-loop column; any other column in that same row has at least two
incoming neighbors. Hence all one-neighbor columns have distinct such
neighbors, and the state transversal matroid is simple. Its dual has no
one- or two-element cocircuit. The full state set minus any one location
is robust, since every source cycle has at least four vertices. Every
actuator location is individually avoidable.

## Both solution maps

An original robust set satisfying the group quotas lifts by selecting
its original representatives. Row coextensions preserve its ranks;
parallel extensions preserve the ranks of the retained representatives.
All source quotas are already satisfied. Cost is unchanged.

Given an arbitrary robust target placement, normalize each source cycle.
Keep its selected a_e, remove all selected fillers and other new source
loops, and add a_e until two group representatives are selected. At
least two source states were selected, so cost does not increase.
Removing loops and adding elements preserves robust matching.

Project selected parallel copies to their support in M'. If some support
element is a coloop of that restriction, a robust parallel selection
contained at least two of its copies. The coloop-protection lemma adds
at most one original ground element per such coloop, without removing
any support element, and returns a no-larger robust support. All original
group witnesses are retained. Finally reverse the row splits. Their
exchanges also retain selected original elements. Thus no selected
original representative is lost, all quota choices survive, and the
result is an original feasible solution of no greater cost. This proves
the optimum identity, including when the matroid has rank zero.

## Size and girth

Use q for the incidence count after compression, h<=n for the grouped
element count, and b for the number of incoming-column splits. Each old
column has degree d_v>=1 and receives at most one source attachment, so

    a <= q-k,       b <= sum_v(d_v-1)=q-n.

Before leaf completion the core has

    N_0 = 2h+n+3k-2+2a+2b <= n+k+4q-2 <= 6q,
    E_0 <= 3*N_0/2 <= 9q.

For any positive odd L, pad this core with the source-cycle port-moving,
private-gate and one-way-arc operations in `source-scc-cover-hardness.md`.
Every bidirected core edge is a private gate, as those operations require.
Up to relocating source ports, the underlying graph is its uniform
L-subdivision. New source states are loops; other additions are parallel
copies or loops. The preceding reverse map is unchanged.

Then complete each old leaf with a bidirected cycle of length 4L instead
of four, attached by one bridge. These cycles give the same parallel/loop
extension, and every cycle in the whole graph has length at least 4L.
There are at most n old leaves, giving

    N = N_0+(L-1)E_0+4L*leaves <= 13qL.

For fixed delta choose integer a>=max(1,(1-delta)/delta) and the smallest
power of three L>=(13q)^a. The output size is polynomial and
N^(1-delta)<=L<=4L. Degrees, bipartiteness and individual avoidability
persist.

## Certificates and scope

`quota_cotransversal.py` constructs the graph and both solution maps;
`sparse_control_graph` accepts an arbitrary feasible control source.
`verify_quota_cotransversal.py` checks redundant and full-rank presentations,
partial and multiple quota groups, every feasible directed source through
three states, rank-zero sources, large-row gates, optional leaf copies,
arbitrary target placements and odd lengths three and nine. Direct
state-plus-input matchings check the selected-actuator-failure model.
Full original-representative rank tables are compared with the source
dual ranks. Other checks cover the parallel rank projection, quota and
original-choice preservation, exact smallest target optima, all source
SCCs, connectivity, degree/parity, simplicity, individual avoidability and
exact girth/size formulas. The certified nontransversal rank-five example
is included with a proper quota subset, so the source need not be an
assignment matroid.
The current certificate set has 618 instances, 68,600 candidate target
placements, 2,001 forward maps, 47,491 reverse maps, 12,326 exact rank
comparisons and 471,271 direct input-failure matching tests. These include
all 220 feasible directed graphs through three states and two larger
positive-deficiency sources with multiple source SCCs and ungrouped
states. The largest padded fixture has 860 states. On targets through
40 states the input-first augmented matching verifier is also compared
with the existing independent state-first augmented matcher.

This composes the earlier row-splitting and parallel-completion exchanges
with source SCCs and a connector of dual loops. Transversal/strict-gammoid
duality is classical; see [Fink--Olarte, Section 7](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/jlms.12505).
The prior [robust-assignment model](https://arxiv.org/abs/1703.06074) and
[fault-tolerant matroid bases](https://doi.org/10.4230/LIPIcs.ESA.2025.83)
do not themselves state this control realization. The sources reviewed
have not established priority for the combined theorem. Its consequence
is structural: any future stronger general hardness ratio can be carried
to this sparse class without a cost baseline or approximation loss.
