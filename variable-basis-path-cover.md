# A variable-basis reduction through path cover

This direction uses the same unit-cost dedicated one-failure model, with
distinct actuator states. It escapes the forced-basis normalization of the
cover construction, but does **not** improve the general 469/468 lower bound.
It establishes an exact source/target equivalence and identifies a further
obstruction to reaching 3/2 through this particular family.

Let H=(W,F) be a simple undirected graph, |W|=n, with no isolated vertices.
Write pc(H) for the minimum number of vertex-disjoint paths partitioning W;
singleton paths are allowed. Allowing cycle components as well gives the
same optimum, since a cycle can be opened into a path without changing the
number of components.

## State construction and exact optimum

Make a state b_v for each v, and two states t_e,c_e for each edge e=uv.
Add the bidirected edges

    c_e--b_u,    c_e--b_v,    c_e--t_e.

There are n+2|F| states and no state self-loops. The graph is bipartite and
is strongly connected when H is connected. Its maximum degree is
max(3,Delta(H)). The matching deficiency is n.

**Theorem.** The robust actuator optimum of this state graph is

    OPT = n + pc(H).

Both directions below give polynomial solution maps: any robust placement
of size s yields at most s-n paths/cycles in H, and any path/cycle cover with
p components yields a robust placement of size at most n+p.

The target matching transversal matroid has one row per c_e, adjacent to
{t_e,b_u,b_v}, plus free right-center columns. Its rank on the target states
is |F|, using c_e->t_e; transposing this matching independently covers all
right centers. Thus center states are loops in the dual matching matroid.

On target states that dual matroid has the following assignment
representation: one job per source vertex v, a private resource b_v for
that job, and a resource t_e able to serve either endpoint of e. To verify
duality directly, a candidate n-element basis A has its selected edge
resources matchable to the jobs whose private resources are absent iff
the complementary targets can be matched to the c_e rows: every absent
t_e forces c_e->t_e, leaving exactly the selected-edge rows for the absent
private vertices. This is a basis-by-basis equivalence.

This assignment matroid is the classical bicircular matroid of H with one
additional graph loop at every vertex: edge subsets are independent iff
every component has at most one cycle. Here these graph loops represent
the private **actuator states** b_v, not self-loops of the state digraph.
For standard bicircular/transversal background see [Bonin, *An Introduction
to Transversal Matroids*, Section 6](https://bpb-us-e1.wpmucdn.com/blogs.gwu.edu/dist/3/152/files/2016/04/TransversalNotes-2ltphpb.pdf).

For a selected resource set S, form the spanning graph K consisting of its
selected ordinary edges and selected private loops. A set spans the
bicircular matroid precisely when each component of K is cyclic. It survives
any one resource deletion iff

1. every vertex has degree at least two, counting a loop twice; and
2. no component is a cycle, including a singleton loop.

Necessity follows by deleting a degree-one edge, or any edge of a cyclic
component with exactly as many edges as vertices. For sufficiency, a
connected graph of minimum degree two that is not a cycle has at least two
independent cycles. Deleting a nonbridge leaves it connected and cyclic.
If a bridge is deleted, neither side can be a tree: such a tree would have
at least two leaves, but only the bridge endpoint can have lost a unit of
degree; a singleton tree side would already have had degree one before the
deletion. Deleting a loop is a nonbridge deletion. Thus every resulting
component remains cyclic.

Each component of H has positive matching rank, so robust spanning already
supplies at least two actuators in its corresponding state component.
Accessibility adds no cost, including when H is disconnected.

## Polynomial extraction: the lower bound

Drop selected center states, which do not help matching. Let s be the
remaining resource count, and F_S the selected ordinary edges. Build a
vertex-disjoint path/cycle cover in F_S, initially all singletons, and apply
these moves until none applies:

- Join endpoint-to-endpoint paths across an edge.
- Join a path endpoint to a cycle, opening the cycle at that vertex.
- Join two cycles across an edge, opening each at its incident vertex.
- Close a path of at least three vertices when its endpoints are adjacent.

The first three moves reduce the number of components. The last preserves
that number and reduces the number of path components. This lexicographic
measure bounds the number of moves polynomially; no optimum path-cover
oracle is used.

Let the resulting cover have p path components and q cycle components.
It uses n-p ordinary edges, so the number of selected resources outside
its edges is s-n+p.

Each nontrivial path needs an outside selected resource at each endpoint,
by the minimum-degree condition. These two resources are distinct, since
an ordinary edge between its two endpoints would permit closing it. A
singleton also has at least two distinct outside resources: at most one is
its private loop, and that loop alone would be a forbidden cycle component.
No resource can serve endpoint charges of different paths, since their
endpoints could then be joined.

Every cycle needs an incident outside resource, since it is not a whole
component of K. Such a resource cannot simultaneously serve a path endpoint
or a different cycle, since a corresponding join move would apply. It may
be a private loop, an internal chord, or an edge to a path's internal vertex.
Thus all charges are distinct, and

    s-n+p >= 2p+q,    hence    p+q <= s-n.

This proves the lower bound and the polynomial extraction map. Opening the
q cycles gives a path partition with the same number of components.

## Lifting a cover: the upper bound

First apply the same improvement moves to the supplied cover, using H's
edges; this never increases the number of components. For each nontrivial
path select its edges and the private loops at its endpoints. For each
cycle select its edges and one private loop on it. These choices cost
|V(component)|+1 each and satisfy the robust graph conditions.

The remaining singleton vertices are independent in H; otherwise their
paths could be joined. Every singleton v therefore has a neighbor u in a
nontrivial path/cycle component. Add edge vu and the private loop at v,
costing two, equal to |{v}|+1. This adds a cyclic leaf side across a bridge
to a component that is already robust, and preserves both graph conditions.
The total cost is at most n+p. Applying this to an optimum path partition
and combining with the lower bound proves the exact formula.

## No forced actuator, and failure of the private-basis normalization

If H has minimum degree at least two, **every individual state can be
omitted by a robust placement**. The target transversal matroid is simple:
different t_e have different sole neighboring rows, and each b_v has at
least two neighboring rows, so every pair of distinct targets is matchable.
Free center columns add no circuits. Thus the dual has no two-element
cocircuit and no mandatory element, by the criterion in
`hardness-obstructions.md`. Strong connectivity makes accessibility
automatic in this claim when H is connected.

This is more than replacing forced states by parallel copies. The evident
private basis B={b_v} cannot be included without loss: if B must be present,
the additional edge resources must form an edge cover of H. On an n-cycle,
this costs n+ceil(n/2), whereas unrestricted OPT=n+1. Consequently the
cost-preserving normalization used for the old gadgets is unavailable for
this basis. Choosing a suitable basis is part of the optimization.

## Parallel copies and degree-three realization

This family nevertheless has a different, variable-basis normalization:
arbitrary parallel copies of **any** resource leave OPT=n+pc(H) unchanged.
Parallel copies here mean copies in the matching matroid: they share one
rank contribution, even if several selected copies protect it from failure.
They are not independent assignment resources with the same endpoints.

Project a robust placement to its original support P. Each coloop of M|P
must have at least two copies. The general exchange lemma in
`hardness-obstructions.md` adds at most one original resource per such
coloop, producing an original robust set of no larger cost. This applies
because every job has an edge option, so the original matroid has no
coloop. It does not normalize to the private basis B. Applying the path/cycle
extraction to this completed support gives the polynomial reverse map.
`cover_from_parallel` implements precisely these two steps.

Now split a high-degree target b_v into a path of target copies, separated
by new degree-two centers, attaching just one old incidence to each target
copy. The old centers c_e retain degree three. This operation makes every
state's degree at most three and preserves the optimum. The relevant
presentation identity is slightly more general than edge subdivision:
redistribute any subset of rows adjacent to a target v to a new target
alpha, and add a row adjacent to {v,alpha}. Contracting alpha recovers the
old transversal matroid; deleting both targets has the old rank after
deleting v. Hence the dual extension is precisely a parallel copy of v.
The matching proof is the same as for the single-incidence operation in
`hardness-degree-threshold.md`, moving whichever redistributed row matched
alpha back to v if the new row was matched to v.

`sparse_realize` implements this split. For m=|F| its graph has 6m-n states
and 7m-2n edges. Strong connectivity and bipartiteness are preserved. The
parallel stability argument proves OPT remains n+pc(H), with no forced
actuators when H has minimum degree two. Subdividing its center--target
edges by odd paths only adds further parallel copies and matroid loops.
The same polynomial padding calculation as in the degree-threshold note
therefore achieves girth at least N^(1-delta), for every fixed 0<delta<1.
Both directions still have polynomial solution maps: project a placement,
extract its cover by the charging argument, then lift that cover.

## What hardness transfers, and what this family cannot prove

Make a (1,2)-TSP instance on H's vertex set by assigning cost one to edges
of H and two to nonedges. Let its optimum be C. Concatenating a path cover
gives a tour of cost at most n+pc(H). Conversely, deleting the cost-two edges
of any tour leaves that many paths; if there are no such edges it leaves a
Hamiltonian cycle. Therefore

    C <= n+pc(H) <= C+1,

with equality on the left unless H has a Hamiltonian cycle. In particular
the exact source objective here is not the number of cover paths by itself;
the unavoidable rank cost n must be retained.

To obtain a hardness reduction with **no mandatory actuator**, start from
an arbitrary (1,2)-TSP source H and add two adjacent universal vertices to
form J. It is connected and has minimum degree at least two. Write C for
the original TSP optimum and R for the robust optimum constructed from J.
Then

    C <= R <= C+3.

For the upper bound, remove cost-two edges of a source optimum tour, or open
its Hamiltonian cycle. The resulting at most max(1,C-n) paths extend to J
by appending the universal vertices to a path. Lift that cover to actuators.
For the algorithmic lower direction, extract a path/cycle cover of J from
any robust placement S, concatenate it into a (1,2)-tour of cost at most
|S|, and shortcut the two universal vertices. All (1,2) costs are metric,
so shortcutting cannot increase the tour cost.

The published [Karpinski–Schmied 535/534 lower bound for (1,2)-TSP](https://arxiv.org/abs/1201.5821)
therefore transfers as 535/534-epsilon hardness for this unforced family:
the additive three is absorbed for sufficiently large n, and bounded n
can be solved by constant-size exhaustive search. This is weaker numerically
than 469/468, and is a different reduction rather than an improvement of
that bound. The source constant is also restated in [Behnezhad, Roghani,
Rubinstein and Saberi,
*Sublinear Algorithms for TSP via Path Covers*, ICALP 2024, Section 1](https://drops.dagstuhl.de/storage/00lipics/lipics-vol297-icalp2024/LIPIcs.ICALP.2024.19/LIPIcs.ICALP.2024.19.pdf).

It also gives a precise obstruction: the known [Berman–Karpinski
8/7 approximation for (1,2)-TSP](https://doi.org/10.1145/1109557.1109627)
yields at most (8/7)OPT+1 actuators in this family by deleting the returned
tour's cost-two edges and lifting its cover. Bounded-size exact search
absorbs the additive one, giving an 8/7+epsilon algorithm for each fixed
epsilon>0 when the source construction is supplied. This already rules out
a larger gap reduction through this source objective, without requiring
recognition of an unlabeled subdivided target. Thus this family cannot support a
3/2 inapproximability threshold unless P=NP. To pursue that target, a next
candidate must use larger resource neighborhoods or a different gammoid,
rather than merely a more complicated source graph H.

## Verification and priority

`variable_basis.py` implements both polynomial maps.
`python verify_variable_basis.py` checks every isolate-free labeled source
graph through five vertices (814 graphs), every robust resource subset on
each, the exact path-partition optimum, both solution maps and absence of
mandatory resources when minimum degree is two. For sources through three
vertices, explicit state-plus-input matching checks the matroid realization
and each actuator failure. Additional fixed-seed six-vertex cases are checked.
Spanning supports with their coloops duplicated are checked as well;
600 full-state rank queries verify the sparse and subdivided realizations
and their reverse solution maps. General private-basis neighborhoods and
the associated hyperforest correction are in `private-basis-hyperforest.md`.
The bicircular definition and the TSP/path-cover relationship are classical;
no literature-priority claim is made for their combination here.
