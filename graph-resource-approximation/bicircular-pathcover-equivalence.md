# Robust graph-resource assignment and path covers

## Scope and statement

Let G be a finite simple undirected graph with no isolated vertices, n=|V(G)|. Every graph edge is a resource that can serve either endpoint job. Add one distinct singleton resource at each vertex; equivalently, add one loop at every vertex. All resources cost one. We must choose resources so that every job can be matched after deletion of any one selected resource. There are no mandatory selected resources or additional quotas in this statement.

Write pc(G) for the minimum number of vertex-disjoint paths covering all vertices, allowing singleton paths. Then

    OPT_robust(G plus singleton resources) = n + pc(G).

Both inequalities have polynomial constructive solution maps. The same identity holds without singleton resources when G has minimum degree at least two and no connected component that is just a cycle.

This is a subclass result, not an algorithm for arbitrary resource neighborhoods or arbitrary source quotas. Transferring it to control graphs requires the already established assignment completion and realization maps, with any accessibility quotas redundant. No claim about general 3/2 optimality follows.

## Graph criterion for robust assignment

The incidence transversal matroid is the bicircular matroid. Its rank on a resource set equals the number of vertices minus the number of tree components of the selected spanning multigraph. Consequently selected resources are robust precisely when every selected component has minimum ordinary multigraph degree at least two and has more edges than vertices. Loops contribute two to ordinary degree.

Proof of the converse direction in this criterion: deleting a nonbridge keeps its component connected and leaves at least one cycle, since the initial component has at least two independent cycles. Deleting a bridge cannot produce a tree component: on a side with v vertices, the original minimum-degree condition gives internal degree sum at least 2v-1, exceeding the 2v-2 of a tree. Thus each resulting component retains a cycle and the matching rank remains n.

For necessity, a vertex with degree zero or one either cannot be matched or is left isolated by deletion of its sole incident edge. A unicyclic component becomes a tree after deletion of a cycle edge. (A single loop at an isolated vertex is such a unicyclic component.)

## Reverse map: at most one path per unit of excess

**Lemma.** A connected multigraph H with minimum degree at least two and |E(H)|>|V(H)| has a vertex-disjoint path cover with at most |E(H)|-|V(H)| paths. Loops and parallel edges are allowed, but paths are ordinary simple vertex paths, including singleton paths.

Start with a minimal connected bicyclic subgraph K of H. Such a subgraph is a theta, a figure-eight, or two cycles joined by a path (a dumbbell), with loop and parallel-edge degeneracies allowed. Each has a Hamilton path through its own vertices. Thus its excess is one and one path covers its vertices.

Suppose K does not yet cover all vertices. There is a component C of H-V(K) adjacent to K.

- If C is acyclic, choose a path through C with both end attachments in K. Such a path exists: every leaf of the tree C has an attachment to K, because H has minimum degree at least two. A singleton C has two attachments. The added structure is an ear with at least one new vertex.
- If C contains a cycle, choose an attachment from K and a shortest path through C to a cycle, and add that path together with the cycle. This is an attached lollipop. A loop is permitted as the terminal cycle.

Each addition has one more new edge than new vertex and its new vertices alone admit a Hamilton path. Append that path as a new member of the vertex-disjoint cover, retaining all previous cover paths. The maintained subgraph remains connected and has minimum degree at least two. Iterate until it spans H.

If p pieces were used, the final subgraph has excess p. Since it is a spanning subgraph of H, p <= |E(H)|-|V(H)|. The construction is polynomial by ordinary cycle/path searches. Apply it independently to each robust selected component and sum to obtain

    pc(G) <= |S|-n.

Deleting loops from the produced vertex paths changes no paths, so these are paths in the original simple graph G.

## Forward map for singleton completion

Start with any vertex-disjoint path cover with k paths. Repeatedly merge two distinct paths when one endpoint of each is adjacent; this does not increase k. A singleton vertex counts as an endpoint. Thus after merging, endpoints of distinct paths are nonadjacent.

For each nontrivial path, select its edges and the two singleton resources at its endpoints. This is a dumbbell with two loops, hence robust, and has one more resource than vertices.

For each singleton path at v, select its singleton resource and one incident graph edge vw. Such an edge exists because G has no isolated vertices. Its other endpoint w is internal to a nontrivial cover path, because all endpoints of distinct cover paths are nonadjacent. This attaches v and its loop to an already robust component. The selected edges for different singleton paths are distinct, and none was a cover-path edge.

The number of resources is (n-k)+2k=n+k. Every resulting component is robust by the graph criterion. Therefore OPT <= n+pc(G), proving equality with the reverse map.

## Forward map without singleton completion

Assume minimum degree at least two, with no cycle connected component. Again merge adjacent endpoints of different paths. Add one unused graph edge at each endpoint; at a singleton vertex add two distinct incident edges. The only possible duplicate added edge closes a nontrivial original cover path into a cycle. If c such duplicates occur, the resulting 2-edge-cover has exactly n+k-c edges.

Every cycle component of this 2-edge-cover is a single original cover path closed by a duplicate edge. Indeed, an added edge entering an internal path vertex would raise its degree to at least three; an added edge between endpoints of different paths is excluded by the endpoint-merging normalization. Thus there are at most c cycle components.

For every remaining cycle component, add a chord or an edge to a different current component. One exists because the corresponding component of the original graph is not a cycle. An added chord makes that component bicyclic; joining two cycles or joining a cycle to a multicyclic component also removes all affected cycle components. At most c added edges are needed. The final graph satisfies the robust criterion and has at most n+k edges.

## Exact relation to (1,2)-TSP

Give pairs that are edges of G distance one and all other pairs distance two. Let TSP12(G) be the minimum Hamilton-tour cost on this complete metric graph (n>=3).

If G is Hamiltonian, TSP12(G)=n and pc(G)=1. Otherwise

    TSP12(G) = n + pc(G).

To see this, a tour with a>0 distance-two edges splits into a paths when those edges are removed, hence pc(G)<=a and its cost n+a is at least n+pc(G). Conversely, joining a minimum path cover cyclically yields a tour costing at most n+pc(G); for one path the closing edge has cost at most two.

Thus robust optimum equals TSP12(G), except that it is one larger when G is Hamiltonian. The maps are polynomial on supplied solutions:

- Robust solution -> the excess path cover above -> a tour of cost at most the robust solution size.
- Tour with at least one distance-two edge -> its path cover -> a robust solution of size at most the tour cost.
- An all-distance-one tour explicitly certifies a Hamilton cycle and gives an optimum robust solution of size n+1.

Therefore an alpha approximation for (1,2)-TSP directly yields an alpha approximation for this robust subclass, with no ratio loss. Conversely a robust approximation yields a tour approximation up to the additive-one optimum discrepancy, which disappears on large instances. This connection does not turn general metric-TSP hardness into (1,2)-TSP hardness.

## Verification and literature caveat

The independent enumerator audit_bicircular_pathcover.py compares exact path-cover DP against exact robust subset enumeration on singleton-completed sources: all no-isolate labelled simple graphs through five vertices and 100 seeded graphs on six vertices. It checked 901 feasible sources. A separate run checked 423 eligible uncompleted sources through seven vertices. These finite checks supplement the constructive proof.

The primary algorithmic source is Berman and Karpinski's 2006 extended-version theorem, cited with its precise path-edge objective below. The stronger bound is a theorem transfer; this directory does not implement their approximation algorithm. No current-best or priority claim is made.

**Source distinction:** [Karpinski, Lampis and Schmied](https://arxiv.org/abs/1303.6437) give the 123/122 lower bound for general metric TSP, while their introduction separately reports 535/534 for (1,2)-TSP. These are different input classes. Consequently this equivalence alone does not prove hardness stronger than the existing 843/838 actuator-placement bound.

## Two mandatory singleton resources

Suppose additionally that the singleton resources at two distinct vertices s,t must be selected; all other singleton resources remain available and optional. Let pc_{s,t}(G) be the minimum path-cover cardinality subject to both s and t being endpoints, possibly of the same path; singleton vertices count as endpoints. Then

    OPT_mandatory(s,t) = n + pc_{s,t}(G).

For the reverse map, choose the initial bicycle in each robust selected component as follows. If it contains both mandatory loops, take those two loops and a simple path joining their vertices. The resulting initial Hamilton path has endpoints s,t. If it contains one mandatory loop, take that loop, any other cycle, and a shortest path connecting them (possibly a figure-eight sharing the loop vertex); its Hamilton path has the mandatory vertex as an endpoint. If it contains neither, use any bicycle. Append ears/lollipops exactly as before; these add disjoint paths and never change endpoints of already constructed paths. The same excess bound gives an anchored path cover with at most |S|-n paths.

The forward map works without merging endpoints. Give each nontrivial cover path its two endpoint singleton resources. Give each singleton cover vertex its own loop and any incident nonloop edge. Every component containing a nontrivial path is already multicyclic; any remaining component has at least two singleton vertices and a loop at each, so is also multicyclic. No vertex is a tree leaf. Duplicated choices only reduce cost below n+k. Both mandatory resources are selected because the required vertices are endpoints. This proves the upper bound, hence equality.

This proof applies to two mandatory singleton resources with all other singleton resources optional. It is not a statement about arbitrary forbidden singleton locations or arbitrary quotas. The independently established 4/3 method may have broader such scope.

The anchored extension is additionally checked by audit_anchored_pathcover.py on 499 graph/anchor-pair instances (all no-isolate labelled sources through four vertices, plus seeded five-vertex sources), comparing an independent permutation-based anchored-path DP with exhaustive robust selected-resource search. Each audit has a main guard; importing its helpers performs no enumeration.

## Exact 8/7 approximation with two mandatory singleton resources

**Corollary.** For a simple graph G without isolated vertices, with one available singleton resource per vertex and mandatory singleton resources at distinct s,t, there is a polynomial-time 8/7 approximation.

### Primary algorithmic source

Berman and Karpinski, *8/7-Approximation Algorithm for (1,2)-TSP (Extended Version)*, ECCC TR05-069 revision 2 (June 2006), Theorem 2, supplies a polynomial-time approximation for maximum path-cover **edge count** L satisfying L_alg >= (6/7)L*. The paper defines the objective as the number of edges in vertex-disjoint simple paths, with isolated vertices counted as zero-edge paths. Its Section 5 explicitly states the equivalent inequality q <= N/7 + 6p/7 for the returned versus optimal path counts on N vertices. There is no bounded-degree, connectedness, or minimum-degree restriction on this theorem.

Primary source:
https://eccc.weizmann.ac.il/report/2005/069/revision/2/download
Theorem 2 is on printed page 2; the objective and inequality are checked in Section 5, printed pages 18–19. This is a cited algorithmic theorem, not a claim to have independently re-proved or implemented its full local-search algorithm.

### Pendant-vertex reduction

Add two new degree-one vertices s',t' to G, with edges ss' and tt', and call the resulting graph G'. Put N=n+2 and p=pc_{s,t}(G).

Every anchored path cover in G extends to a path cover in G' of the same cardinality: append the new leaf at each required endpoint. If an anchor is a singleton, appending its leaf simply creates a two-vertex path. Thus pc(G')<=p.

Conversely, normalize any path cover of G' to contain both pendant edges. If ss' is missing, s' is an isolated path. If s has cover degree at most one, add ss'; otherwise delete one of the two cover edges incident to s and then add ss'. The result remains a collection of paths and its edge count does not decrease. Repeat at t. Since s and t are distinct, normalizing at t cannot delete ss'. Both pendant edges are now selected. Remove s',t' and their incident edges. This preserves the number of paths because each removed leaf lies on a path containing its original anchor; no entire path disappears. It also leaves both s,t with cover degree at most one, so they are endpoints. Hence p<=pc(G') and equality follows.

### Approximation calculation

Run the Berman-Karpinski maximum-edge path-cover algorithm on G'. Its optimal edge count is N-p. The returned edge count L satisfies

    L >= 6(N-p)/7.

Normalize pendant edges and remove the two new leaves as above. Since normalization never decreases edge count, the recovered anchored cover has q paths with

    q <= N-L <= N/7 + 6p/7.

Map this anchored cover to a robust resource set by the constructive theorem. Its cost is at most

    n+q <= (8n+2+6p)/7 <= 8(n+p)/7,

where the last inequality is exactly 2<=2p, true since the original graph is nonempty and p>=1. Finally n+p is the optimum mandatory-resource cost. This proves the exact 8/7 guarantee.

The two added leaves exist only in the algorithm's auxiliary path-cover graph. They do not add selected resources, actuator locations, costs, mandatory locations, or source-SCC quotas to the original instance. Other singleton resources must remain available for the constructive forward map; arbitrary forbidden resources and arbitrary quotas are outside this corollary.


## Exact transfer to the four-source-state control construction

Consider the state graph constructed from a simple graph H with n jobs, ordinary resources for its edges, physical singleton resources at distinct s,t, completion singleton states at all jobs, and source states A,B,F_A,F_B. Its assignment interpretation maps A and B to additional singleton copies at s and t. The fillers F_A,F_B and all center states are loops of the dual matching matroid. Accessibility requires at least two of the four source states. The construction and independent state-graph verifier are described in [the 4/3 note](four-thirds-graph-resources.md), Section 5.

Let Y be the set of center and filler states, all loops of the dual matching matroid. Any robust actuator placement S can be normalized to

    S' = (S minus Y) union {A,B}.

The source quota gives |S intersect {F_A,F_B}| >= 2 - |S intersect {A,B}|. Thus the number of added source states is at most the number of removed fillers, and |S'| <= |S|. This normalized placement selects A and B and no filler or center states. If A or B is missing, the source quota supplies a selected filler that can be replaced by that missing source state. Replacing a matroid loop by an arbitrary new element preserves robust spanning: after failure of the new element, the old set with its loop removed still spans; after failure of any other element, its old surviving spanning set remains present except for a useless loop. Each replacement preserves the source quota. Once A and B are selected, delete all remaining fillers and centers. These deletions preserve robust spanning and accessibility.

The normalized selected resources contain mandatory singleton copies A and B, together with any selected ordinary resources, completion singletons, and extra physical singleton copies at the anchors. The reverse anchored bicycle-and-ear construction allows loops and parallel resources. Start with the two specified loops if they lie in the same component, or with the specified loop and another cycle in each relevant component. In the parallel-loop degeneracy the initial path may be the singleton anchor. Thus the extra copies do not invalidate the extraction, and the resulting anchored path cover satisfies

    pc_{s,t}(H) <= |S| - n.

Conversely, the forward anchored construction uses A and B as the selected endpoint singleton resources at s,t, and completion singletons at any other endpoints. Its assignment is robust and its source quota is satisfied. Therefore the optimum over all allowed actuator locations, including every auxiliary state, is exactly

    OPT_control = n + pc_{s,t}(H).

The pendant reduction and Berman-Karpinski theorem consequently give an exact 8/7 approximation for this entire control family. It remains a subclass result and does not improve the general arbitrary-state-graph ratio.

## Pendant-reduction audit

`audit_pendant_anchors.py` exhausts 259 base-graph/anchor instances through four vertices. It checks 10,220 augmented path covers, confirms the two-edge optimum shift, verifies both pendant-edge normalizations and preservation of path count on deleting the leaves, and checks the final ratio for all 703 covers meeting the cited 6/7 path-edge guarantee. This verifies the reduction independently of an implementation of the cited approximation algorithm; no such implementation is included.
