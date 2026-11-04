# A degree threshold for single-actuator robustness

Model throughout: unit-cost dedicated actuators, at most one per state, any
one selected actuator can fail, and independently free structural parameters.
Outdegree counts distinct outgoing state neighbors, including a self-loop.
Girth and planarity below refer to the underlying **undirected** graph; a
bidirected graph always has directed two-cycles.

## Results

1. On arbitrary directed graphs of maximum outdegree two, 1-RMSCP-d is
   polynomial-time solvable. An exact formula is given below, including source
   SCC accessibility and infeasibility.
2. It is APX-hard on strongly connected, loopless, bidirected, bipartite graphs
   of maximum undirected degree three. On this class, approximation within
   `469/468 - epsilon` is NP-hard for every fixed positive epsilon for which
   the factor exceeds one, unless P=NP.
3. The hardness in item 2 persists at every prescribed fixed undirected girth.
   With planarity also required, the decision problem remains NP-complete;
   no APX-hardness assertion is made for the planar subclass.

Thus two versus three is an exact complexity threshold for a bound on the
maximum state outdegree. The constant is a modest improvement over 694/693,
not a 3/2 inapproximability result. The principal result is the structural
threshold and the cost-free port/subdivision construction. Literature priority
of these structural statements remains unconfirmed after a targeted review.

The starting reference is [Poças, Pequito and Ramos, Automatica 190 (2026),
113074](https://doi.org/10.1016/j.automatica.2026.113074), Theorems 3–4.
Their Set Cover construction charges two forced actuators per universe
element: OPT=2|U|+tau. The construction here retains that cost relation and
replaces high-degree ports without increasing the actuator baseline. For the
source approximation gap, [Chlebík and Chlebíková, *Complexity of approximating
bounded variants of optimization problems*, TCS 354 (2006), 320–338](https://doi.org/10.1016/j.tcs.2005.11.029)
give 53/52 hardness for Vertex Cover on 4-regular graphs. A
[public author manuscript](https://pure.port.ac.uk/ws/portalfiles/portal/1887750/3DM_JOURNAL_revision_old.pdf)
states this on p.3 and in Corollary 18. This use of the stronger source gap
also improves the numerical constant in the original unbounded-degree
construction; it is not by itself a new gadget principle.

## Outdegree at most two

Let T be the state transversal matroid on right states, M=T*, and C_1,...,C_l
the source SCCs. A placement S is robust precisely when every S-s spans M
(also S itself spans) and |S intersect C_j|>=2. This is the characterization
in `bounded-core-theorem.md`.

Make a multigraph H on the state vertices. A row of the state adjacency
matrix with two neighbors u,v supplies an edge uv; a row with sole neighbor
v supplies a loop at v; an empty row supplies nothing. Distinct rows supply
distinct edges, possibly parallel. The rows, rather than the states, label
these edges. Matching a subset of right states to distinct left rows is
exactly assigning these vertices to distinct incident edges of H.

For a component Q of H:

- If Q is a tree, any proper subset of its vertices can be assigned distinct
  incident edges: root at an omitted vertex and use each selected vertex's
  parent edge. All vertices cannot be assigned, since |E(Q)|=|Q|-1.
- If Q contains a cycle, take a spanning connected unicyclic subgraph.
  Assign the cycle vertices their successive cycle edges, and each remaining
  vertex its parent edge toward the cycle. This also works for a loop and a
  two-edge cycle of parallel edges. Thus every vertex can be assigned.

The edge sets for distinct components are disjoint. Consequently T is a
direct sum of U_(|Q|-1,|Q|) for tree components Q and free matroids on the
other components. Its dual is a direct sum of U_(1,|Q|) on the tree components
and loops elsewhere, so the matching deficiency equals their number p.
Let Q_1,...,Q_p be these tree vertex sets. Matching
robustness is therefore **exactly** |S intersect Q_i|>=2 for every i.

A singleton Q_i without a loop means a state has no incoming state edge;
it is a singleton source SCC, so the ordinary source feasibility test detects
it. As in the main repository, the instance is feasible iff every source SCC
has at least two states. On feasible instances all Q_i also have size at
least two.

There are now two disjoint families of quota groups, the Q_i and the C_j.
Form a bipartite multigraph with a node for each group, and an edge labeled
v between Q_i and C_j whenever v belongs to both. Each edge has capacity
one and each group node capacity two. Let h be a maximum b-matching size.
It is computed by integral maximum flow. Then

    OPT = 2p + 2l - h.

For the lower bound, choose two witnesses from each group in any feasible
placement. The witnesses shared between the two families form a b-matching,
so there are at most h shared witnesses. Their union has size at least the
displayed expression and lies in the placement.

For the upper bound, select the states labeling a maximum b-matching. Add
states to complete the Q_i quotas, then the C_j quotas. Availability follows
from group sizes at least two. The initial sum of unmet quotas is
2p+2l-2h, and each addition meets at least one unmet quota. The resulting
size is at most h+(2p+2l-2h); the lower bound gives equality. No repeated
actuator is used. This proves the formula and a polynomial exact algorithm.
It applies with arbitrary source SCCs and self-loops, not just bidirected
graphs. `outdegree_two_solver.py` returns an optimum placement.

## A subcubic reduction with no additional actuator cost

Take a simple Vertex Cover instance H=(W,F), m=|F|>0; discard isolated
vertices. Write tau(H) for its optimum. For every edge e=uv introduce a
center c_e and targets a_e,b_e,x_e, joining c_e to all three targets.

At each source vertex v order its incident edges e_(v,1),...,e_(v,k).
Introduce centers p_(v,i) and targets q_(v,i), with edges

    p_(v,i)--q_(v,i),             1 <= i <= k,
    q_(v,i)--p_(v,i+1),           1 <= i < k,
    p_(v,i)--x_(e_(v,i)),         1 <= i <= k.

Replace each undirected edge by its two directed arcs. Call the center set L
and target set R. Apart from the connectors described next, all edges run
between L and R. Each x_e has degree three; each p has degree at most three;
all other targets have degree at most two, except the degree-one twins.

For each source component choose one vertex v and prefix its port chain by
p_0--q_0--p_(v,1), where p_0 is a center. Join these p_0 centers in a path
across source components. Their degrees remain at most three. These new
center--center edges are bridges, and independently flipping the two colors
of consecutive components makes the whole graph bipartite. It is connected,
hence its bidirected version is strongly connected. The graph has 8m+2c
states for c source components. Prefixes and bridges charge no actuators.

Let B={a_e,b_e:e in F}. Define a state matching on targets by

    c_e -> x_e;    p_(v,i) -> q_(v,i);    p_0 -> q_0.

It saturates L. Independently match all right centers using left a_e, q_(v,i)
and q_0. Thus the matching deficiency is |R|-|L|=2m and B is a basis of M.
More generally the state transversal matroid is T_R direct-sum a free matroid
on the centers: matching target subsets can use only left centers, and any
such matching coexists with the fixed matching of all right centers from
left targets. The connector edges enter only free center columns and do not
change this direct-sum statement. In particular, actuating a center gives
no matching benefit.

**Forced baseline.** Every ordinary feasible placement must actuate at least
one of a_e,b_e, since both have only c_e as an incoming state neighbor.
Single-failure robustness forces both. Hence every robust S contains B,
charging exactly 2m unavoidable actuators.

**Exact remaining constraints.** For e=uv at positions i and j in its
endpoint port orders, put

    R_e = {x_e} union {q_(u,t):t>=i} union {q_(v,t):t>=j}.

Then S is robust iff B is contained in S and S meets every R_e.

For necessity, fail a_e. If S misses R_e, all targets in {a_e} union R_e
are unactuated. Their left neighborhood consists of c_e and the two suffixes
of port centers. It has one fewer vertex than this target set. Hall's
condition fails. The same argument applies on failing b_e.

For sufficiency, the displayed target matching admits an alternating path
from either failed twin to every vertex of R_e: go through c_e,x_e, then
forward along an endpoint port chain if needed. End at an actuated z in R_e
and flip the path. All unactuated targets are matched; right centers retain
their independent fixed matching. A failure outside B leaves the original
spanning basis B. Accessibility holds after every failure because the graph
is strongly connected and at least one of its 2m forced actuators remains.

Select q_(v,k), the terminal of v's chain, for each v in a source vertex
cover. These terminals meet all R_e. Conversely map any selected q_(v,i) to
v, and any selected x_e to either endpoint of e. Ignore centers and prefix
targets, which occur in no R_e. The resulting vertices cover H and their
number is at most |S|-2m. Consequently

    OPT(G(H)) = 2m + tau(H),

with solution maps in both directions and no hidden forced port cost.

For maximum source degree Delta, m<=Delta*tau(H). The reduction is an
L-reduction with constants (2*Delta+1,1). If a rho-approximation returns S,
the extracted source cover has size at most

    rho*tau(H) + 2m*(rho-1)
      <= [1+(2*Delta+1)*(rho-1)]*tau(H).

For Delta=4, a factor below 1+(53/52-1)/9 = **469/468** would contradict
the cited source hardness. APX-hardness also follows from subcubic Vertex
Cover, with constants (7,1).

If H is supplied with a planar embedding, order the ports in its cyclic
incidence order. The port chain is a tree drawable in the vertex disk with
these ports on its boundary. Put the three-leaf element gadget in each edge
disk and use the original incidence corridors. Prefixes and intercomponent
bridges can be drawn in the exterior faces. Thus planarity is preserved.
The NP-completeness of planar Vertex Cover is established in [Garey, Johnson
and Stockmeyer, *Some simplified NP-complete graph problems*, TCS 1 (1976),
237–267](https://doi.org/10.1016/0304-3975(76)90059-1).
Membership in NP follows by checking a state/input matching and accessibility
for each selected actuator failure.

## Preserving the optimum while increasing girth

The following argument is necessary: subdivisions introduce new selectable
states, so preserving the deficiency alone would not preserve the optimum.

First consider a bipartite presentation of T_R with all its left rows
saturable. Replace one presentation edge l--v by l--alpha--beta--v, with
alpha a new target and beta a new center. The new transversal matroid T'
has rank r(T_R)+1. Contracting alpha gives T_R: a matching using l--v can be
replaced by l--alpha,beta--v; one not using that edge can match alpha to
beta. The converse follows by removing alpha from a matching and, if needed,
replacing l--alpha,beta--v by l--v.

Write M'=T'*. Deleting alpha from M' gives M. Moreover deleting both v and
alpha from the presentation lowers its rank by one if v is a nonloop of M,
and by two if v is a loop of M. The dual rank formula therefore gives
r_(M')({v,alpha})=r_M({v}). Their singleton ranks have this same value.
For alpha, this follows also from r(T'-alpha)=r(T_R-v)+1; for v it follows
from M'-alpha=M.
Thus alpha is a **parallel copy of v** in M', or another loop when v is a
loop. The new right center is a matroid loop, as all centers remain freely
matchable. Iterating proves that the subdivided state's dual matching
matroid is a parallel extension of the original, plus loops.

For the reduction above, even unrestricted parallel copies cannot reduce
the optimum. Here is a normalization proof. The original a_e,b_e form a
series pair of M (their pair is a circuit of T_R), so exchanging their names
is a matroid automorphism. Also {a_e,b_e,x_e} is a circuit of M: x_e replaces
either twin in the basis B and no other basis element.

In a parallel extension let A_e,D_e be the selected copies of a_e,b_e.
Their union must have size at least two, since the original two-element
cocircuit expands to all copies of its two members. Normalize one pair at
a time:

- If only one member has copies, replace them by the original a_e,b_e.
  Deleting either new actuator leaves a spanning set: the old support with
  one member spans, and the series-pair automorphism exchanges the members.
- If both have copies and at least one has only one copy, retain just one
  original of each. Failure of the singly represented old member shows that
  the other member alone with the remaining support spans; the automorphism
  proves the other deletion case.
- If both have at least two copies, replace their at least four actuators by
  a_e,b_e,x_e. After either twin fails, the surviving twin together with x_e
  spans both twins, by the three-element circuit. After x_e fails, both twins
  remain. Failures elsewhere preserve the original support.

Every step preserves robust spanning and does not increase cost. After all
pairs have been normalized, the full basis B is present. Extra parallel
copies of a state outside B can be discarded: failure outside B leaves B,
and repair of a basis deletion depends only on the support, not multiplicity.
This yields a robust original placement of no larger cost. Original
placements remain robust in the extension. The two optima are equal, and
the normalization is polynomial, providing the required reverse solution
map for approximation hardness as well.

Now replace every center--target edge by a length-three path; leave the
center--center connector bridges alone. Degrees, bipartiteness, connectivity
and planarity are preserved. Every undirected cycle stays inside a source
component, so its length is multiplied by three. Repeat a fixed number of
times to exceed any prescribed fixed girth g. The graph size still grows by
only a constant factor depending on g, and OPT remains 2m+tau(H). This proves
the stated large-girth hardness.

## Exact computational checks

Python 3.10+, standard library only:

    python verify_outdegree_two.py
    python verify_hardness_gadgets.py

The first script exhausts all 15,002 looped directed graphs of maximum
outdegree two on one through four states and checks 500 fixed-seed graphs
on five through nine states, plus four quota-interaction fixtures.
A subset-DP matching oracle checks every actuator
placement, the transversal decomposition and the exact returned optimum.

The second checks all 71 nonempty labeled source graphs on two through four
vertices, every possible one-state repair of each forced twin, solution
lifting/extraction and an exact search for the optimum residual cost. It
also checks K5 and the four- and six-leaf stars, exercising longer ports.
It fully enumerates all actuator placements on the 10- and 18-state targets
of a single edge and a two-edge path. On four source families it enumerates
133,184 matroid support subsets to compare distinct and repeated-actuator
optima, checks 1,025 full-state subdivision rank queries, and checks the
girth multiplication through two rounds. It also exhausts a 15-state
counterexample to the
one-leaf compression in `hardness-obstructions.md`. Tests supplement the
proofs; they do not establish an asymptotic approximation bound.
