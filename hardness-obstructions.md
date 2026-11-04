# Why two natural routes to 3/2 hardness fail

The model is the dedicated, unit-cost, single-selected-actuator-failure
problem. These are obstructions to particular reduction architectures, not
an upper bound for the general problem. The new structural lower bounds
and exact outdegree threshold are in `hardness-degree-threshold.md`.

## Independent flats

Let T be the state transversal matroid and M=T*. For S selected and
I=V\S, matching feasibility says that I is independent in T. Surviving a
failure at s in S says that I+s is independent. Thus matching robustness is
equivalent to I being an **independent flat** of T. Source SCCs additionally
require |S intersect C|>=2. In particular, on a strongly connected graph of
positive deficiency the optimization is exactly the minimum complement of
an independent flat of T; accessibility is automatic for a robust spanning
set, which already has at least two elements.

## Parallel copies do not lower the matroid optimum

Let M be any coloop-free matroid with a polynomial rank oracle. Add arbitrary
parallel copies of its elements and arbitrary loops. The minimum size of a
set spanning after any one deletion is unchanged, with a polynomial reverse
solution map. This statement concerns the matroid condition; additional
accessibility quotas require a separate argument.

Indeed project a robust extended set to its original support P. It spans M.
Every coloop of M|P must have at least two selected copies: deleting its only
copy would otherwise destroy spanning. Let K be this set of coloops, so the
extended cost is at least |P|+|K|.

For b in K that is still a coloop of the current support P', rank(P'-b)=d-1.
Since M has no global coloop, an original e outside P' exists with
rank(P'-b+e)=d. Add e. This protects b; adding elements cannot make any old
noncoloop a coloop, and e is not a coloop because P' already spanned.
At most |K| additions therefore produce a robust original set of size at
most the extended cost. Conversely an original robust set remains robust
in the extension. The rank-oracle procedure is `protect_coloops` in
`hardness_gadgets.py`.

Consequently removing literal mandatory states by parallel duplication does
not by itself escape a forced-basis reduction: the original optimum can
still be recovered. The private-basis example in
`variable-basis-path-cover.md` instead has a genuine loss if that particular
basis is required to be present.

## Removing one forced leaf destroys the Vertex Cover encoding

In the original Set Cover gadget, remove one of the two forced private
leaves for each element. Specialize the source to Vertex Cover on a graph
H=(W,F), with universe F and a set of incident edges for each source vertex.
For each e introduce c_e,a_e,x_e, with c_e joined to a_e,x_e. For each
v introduce s_v,t_v, with s_v joined to t_v and all x_e for e incident to v.
All edges are bidirected. Assume H connected and m=|F|>0.

The tempting claim OPT=m+tau(H) is false. In fact

    OPT = m+1.

The target has matching deficiency m, so every robust placement has size
at least m+1. For the converse choose a spanning tree K of H and put
actuators at every t_v and at a_e for every e outside K. Their number is
|W|+(m-|W|+1)=m+1.

Here are explicit state matching witnesses. With no failure, match the
unactuated a_e for e in K via c_e, and their displaced x_e via distinct
endpoint centers s_v: root the tree anywhere and assign each edge its child
vertex. For every e outside K match x_e via c_e.

If a selected a_f fails, use K+f instead. This is a connected unicyclic
graph, whose edges can be assigned distinct endpoint vertices: orient its
cycle cyclically and its attached trees outward. Match a_e via c_e and x_e
via the assigned s_v for e in K+f, and match other x_e via c_e.

If t_v fails, match it via s_v. Root K at v, and assign each tree edge its
child vertex, obtaining the required distinct centers other than s_v for
the displaced x_e. In all cases right centers are independently matched
using left a_e->c_e and t_v->s_v. These left state vertices remain available
when an actuator fails. The connected state graph supplies accessibility.

For disconnected H without isolated vertices the same proof gives
OPT=m+c(H). Thus the modified gadget cannot encode Vertex Cover even on
connected cubic instances.

Explicit smallest cyclic source: H is a triangle, m=3, tau(H)=2. Number
the three edge gadgets (c_e,a_e,x_e) by (0,1,2),(3,4,5),(6,7,8), and the
three vertex gadgets (s_v,t_v) by (9,10),(11,12),(13,14). Join x_01 to
s_0,s_1, x_02 to s_0,s_2, and x_12 to s_1,s_2, in addition to the private
gadget edges. The placement {1,10,12,14} is robust and optimum, of size 4,
whereas m+tau(H)=5. `verify_hardness_gadgets.py` checks all 32,768 actuator
subsets of this target independently.

## A general obstruction for forced-basis covering reductions

Let M be a coloop-free rank-d matroid, d>0. Consider minimum-cardinality
sets S spanning M after any one deletion. Suppose the set B of **all**
elements mandatory in every such S is itself a basis. This is exactly the
forced-basis architecture of the published cover reduction. Then this
optimization class has a polynomial **5/4 approximation** when its
fundamental circuits are accessible, in particular for state matching
matroids. Consequently a 3/2 hardness reduction confined to this class
cannot succeed unless P=NP.

First, an element e is mandatory iff it lies in a two-element cocircuit.
Since M has no coloops, V-e spans. It fails robustness iff V-{e,f} fails
to span for some f, equivalently a cocircuit lies in {e,f}; no singleton
cocircuit exists, so {e,f} is the cocircuit. Conversely robust spanning
must select both members of a two-element cocircuit. This also shows that
the mate f of a mandatory element is mandatory.

Two-element cocircuits partition B into series classes B_1,...,B_q, each
of size at least two. Indeed, cocircuit elimination makes the relation
transitive, since singleton cocircuits are absent. Hence d=|B|>=2q.

For x outside B, let C(x,B) be its fundamental circuit. Circuit–cocircuit
orthogonality shows that C(x,B) intersect B is a union of whole series
classes: it cannot meet a two-element cocircuit in just one element.
Let F_x be the classes occurring in this union. B+X is robust iff the F_x,
x in X, cover all q classes. To see necessity for b in B, the hyperplane
cl(B-b) contains X precisely when no x exchanges for b. Sufficiency follows
from the replacement basis B-b+x; deletion of an x leaves B.

Write tau for this unweighted cover optimum. The target optimum is exactly
d+tau, but the number of independently encoded requirements is at most d/2.
This explains why merging forced pairs cannot retain arbitrary independent
covering requirements while lowering their baseline below two each.

There is a stronger, algorithmic barrier than this count alone. Make a graph
on the q classes, joining two classes whenever some F_x contains both. Find
a maximum matching of size mu. Choose one covering set for each matched
pair and one for each unmatched class. At most q-mu sets are selected.
Unmatched classes are independent in this graph, so every F_x contains at
most one of them. Therefore tau>=q-2mu and

    ALG <= d+q-mu <= d+(q+tau)/2 <= (5/4)(d+tau),

where the last inequality follows from d>=2q (and tau>=0).
All operations are polynomial. In the state matching setting ranks and
fundamental circuits are found by matching queries; the mandatory set is
also found by testing the polynomially many pairs for cocircuits.

This obstruction assumes that the entire mandatory set is a basis and
that there are no extra source-SCC quotas outside the matroid condition.
It does not cover arbitrary non-forced-basis constructions, cost changes,
or general accessibility interactions. A route to substantially stronger
hardness must escape these assumptions rather than just shorten the
existing private-leaf gadget.

## Reversing the bicircular construction also collapses

There is a different tempting family with a smaller rank baseline. Let H
be a connected simple graph on n vertices and m>n edges. Use its bidirected
incidence graph as the state graph: a center for each vertex and a target
for each edge, adjacent to its two endpoint centers. The state transversal
matroid on targets is the bicircular matroid B(H); the right centers are a
free summand, since a spanning unicycle supplies a matching covering them.
Thus the matching deficiency is m-n. Cubic sources already give degree-three
bipartite, strongly connected targets, without forced actuators.

Nevertheless every such instance has

    OPT = m-n+1.

Take the unactuated edge targets to be any spanning tree of H, together
with all center states. They are independent in the state transversal
matroid, and adding any one omitted target gives a unicycle and remains
independent. This is an independent flat, hence a robust placement. Its
size m-n+1>=2 supplies accessibility. The reverse bound is rank(M)+1,
because deleting any selected element must leave at least rank(M) elements.
Choosing which bicircular side to dualize therefore does not by itself
produce a new hardness source. The private-basis hyperforest objective in
`private-basis-hyperforest.md` is the remaining, more general candidate.
`audit_hardness.py` checks the spanning-tree placement and every actuator
failure on all 388 connected simple source graphs with m>n through five
vertices; the returned full matchings also certify the matching deficiency.
