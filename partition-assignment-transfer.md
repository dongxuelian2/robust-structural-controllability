# Exact transfer with resource-group quotas and genuine choices

The source-SCC method extends beyond a fixed private basis. Let U be n>=1
jobs and E a collection of m unit-cost resources with nonempty job
neighborhoods. Every resource is vulnerable. Suppose the full resource
set is a feasible one-failure robust assignment. Partition E into groups
P_1,...,P_k of size at least two and impose

    |X intersect P_i| >= 2

in addition to robust assignment. Write OPT_P for this constrained
minimum. It is finite because E satisfies both kinds of constraint.

**Theorem.** This problem has a polynomial control-graph realization with
exactly OPT_P actuators, cost-nonincreasing polynomial solution maps in
both directions, and every added state selectable. The graph is loopless
and bipartite, every underlying undirected degree is two or three, both
directed degree bounds are three, and no actuator location is mandatory.
For every fixed 0<delta<1 one may additionally require undirected girth
at least N^(1-delta). If the source job/resource incidence graph is
connected, the target is weakly connected.

There is no additive actuator baseline. A lower bound for this
partition-constrained assignment problem would therefore transfer with
the same approximation factor. The theorem does not itself establish a
new numerical gap or a 3/2 hardness result. Unlike the fixed-basis cover
construction, groups can have three or more functional resources with
different neighborhoods, and a solution map does not fix their choices.

For a concrete separation, use the cycle family in
`variable-basis-path-cover.md`: n vertex jobs, their n private resources,
and n edge resources serving their two endpoints. Put all 2n resources
in one group, whose quota is redundant. The true optimum is n+1, while
requiring all private basis resources costs n+ceil(n/2). The transfer
preserves both facts. Thus no cost-nonincreasing completion to that fixed
basis is available even in this realized family. The cycle family is
easy; this is a separation of architectures, not a hardness claim.

## Source cycles and private routing trees

For each original resource e create a state e and a filler state f_e.
For a group P=(e_1,...,e_s), make the source cycle

    e_1 -> f_(e_1) -> e_2 -> f_(e_2) -> ... -> e_s -> f_(e_s) -> e_1.

Create one physical job sink z_j for each j in U. From the predecessor
f_(e_(i-1)) of resource e_i, attach a private binary routing tree to
precisely its job neighbors. Its root has only one outgoing tree arc.
Use the same internal pair x,y as in `source-scc-cover-hardness.md`:
parent->x, x<->y, and at most two child arcs leaving y. A singleton
neighborhood connects directly to its job sink. Internal trees are
disjoint except for the shared physical job sinks.

Give original resources, physical jobs and internal x states one color;
give all source fillers and internal y states the other. Every arc
crosses colors. This works with one source copy: placing a loop filler
between consecutive resources removes every source-cycle parity issue.

The source SCCs are exactly the displayed cycles. Every other SCC is
a routing pair or a job sink and has incoming arcs. Thus accessibility
is exactly quota two in every expanded source group.

## Matching and the source-filler normalization

In the state transversal presentation, row e contains only f_e, so each
source-filler column is a free summand of T and a loop of M=T*. Likewise
an internal x row contains only y, making y an M-loop. Remove those
free columns with their singleton rows. Permuting the remaining rows
gives the disjoint resource-tree presentation: each original resource
or functional tree state has its own column and its child columns in
its row, and jobs have no outgoing row. Hence the matching deficiency
is n, and the restriction to original resource states is exactly the
source assignment matroid.

Any robust source assignment X lifts by selecting those same resource
states. An assignment matching uses at most one unit from each private
tree, so it lifts to a vertex-disjoint linkage; equivalently its state
matching follows by flipping the tree's private matched pairs. This
works after every original resource failure. The source quotas are
already satisfied and survive each failure by accessibility.

For the converse, first normalize a robust target placement within
each source cycle. Keep every selected original resource. Discard its
selected source fillers, and add original resources of the group until
at least two are selected. At least two cycle states were selected, so
this does not increase cardinality. It removes only M-loops and adds
matroid elements, preserving robust matching. Quotas remain satisfied.

Relax the capacity constraints inside each routing tree. Every remaining
functional tree state then becomes an assignment resource serving its
reachable jobs, a subset of its owning original resource's neighborhood.
Each selected physical job is a singleton assignment resource. A true
linkage yields an assignment, so a robust target placement yields a
robust assignment in this relaxation. The dominated-resource theorem in
`assignment-hardness-transfer.md` removes these added resources without
increasing cost. Crucially that normalization **never removes a selected
original resource**. Therefore it preserves all original group quotas,
including the choices already made within each group. Together with
the lift this proves equality of the optima before degree compression.

The assumption that the original full resource set is robust is needed
by the dominated-resource exchange. It is checkable by polynomially many
matchings; it is not silently inferred from target feasibility.

## Compressing job indegrees while retaining the quotas

Before this step only physical job states can have degree greater than
three. If a job v has more than three incoming state rows, choose two
of them, move their arcs from v to a new x, and add

    x -> y,       y -> x,       y -> v.

The new x has undirected degree three and y degree two; the indegree
of v decreases by one. The y column is a free transversal summand.
After removing it with its singleton x row, the remaining operation is
the standard target split with row {x,v}. Thus x is a parallel copy of
v in M, and y an M-loop. Repeat until every job has degree at most three.
The bipartition and all outgoing degree bounds are preserved.

These parallel copies lie outside every source SCC. After source-filler
normalization, project a robust placement to its support in the old
matching matroid and apply coloop protection from
`hardness-obstructions.md`. It returns a no-larger robust old support
by adding elements, without removing any old support element. In
particular every selected original source resource remains selected,
so all quotas remain true. Now apply the tree-capacity relaxation and
dominated-resource normalization above. This supplies the required
reverse map after degree compression as well.

If q=sum_e |N(e)|, the routing construction has exactly 2q+n states.
A feasible uniform source has every job degree at least two. If job j
has source frequency d_j, compression adds
2*max(d_j-3,0) states. Consequently

    N_0 <= 2q+n+2*sum_j(d_j-2) = 4q-3n < 4q,
    E_0 <= 3*N_0/2 < 6q.

All degrees are at least two. Before compression, columns with just one
incoming neighbor have distinct such neighbors, and all other columns
have at least two. Thus T is simple, or M cosimple. Parallel and loop
extension preserve cosimplicity. Since every source SCC has at least
four vertices, the full state set minus any one state is a robust
placement. Every location is therefore individually avoidable.

## Girth without an actuator surcharge

Apply the three odd-padding operations proved in
`source-scc-cover-hardness.md`: move source-cycle ports while inserting
M-loops, lengthen internal bidirected gate edges, and replace other arcs
by paths of private directed two-cycles. The latter operations add only
parallel copies and M-loops. No original source resource is copied in
its source SCC. Source-filler normalization and coloop protection still
keep the original quota witnesses, after which the same dominated-resource
map applies. The optimum is unchanged.

Up to relocating source ports, the underlying graph is the uniform odd
L-subdivision, with N=N_0+(L-1)E_0<6qL and girth at least 4L. Choose L
as the smallest power of three at least (6q)^a, where fixed integer
a>=max(1,(1-delta)/delta). Then the output is polynomial in q and has
girth at least N^(1-delta). Degrees, bipartiteness and individual
avoidability persist.

## Certificates and scope

`partition_assignment_transfer.py` implements the construction, job
splitting, odd padding and the reverse map. `verify_partition_assignment.py`
enumerates feasible small uniform assignments and every resource partition
whose groups have size at least two. It independently checks original
resource ranks, actual state-plus-input matchings under failures, both
solution directions, source SCCs, simplicity, degree/parity properties,
and complete target optima on the smallest instances. Larger fixtures
exercise arbitrary target selections, job-copy normalization, source
fillers, full-state rank projection and exact girth scaling.
The four-cycle fixture explicitly checks OPT five versus six when its
private basis is required.
The current certificate set has 322 small source/partition instances,
six padded fixtures and that cycle fixture: 24,542 target placements,
416 forward maps, 15,426 reverse maps and 2,632 full-state rank projections.

The matching-only uniform assignment problem and its hardness were studied
in [Adjiashvili--Bindewald--Michaels](https://arxiv.org/abs/1703.06074).
The unit-threshold Activation Edge-Cover barrier for additive fixed-basis
covering is discussed in `hardness-obstructions.md`; it does not apply
to arbitrary partition-constrained robust assignment. No stronger
hardness constant for that source problem is claimed here. Priority of
this exact sparse control realization is unconfirmed beyond the sources
reviewed.
