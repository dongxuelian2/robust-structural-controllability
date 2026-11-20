# NP-completeness at the absolute matching lower bound

## Theorem and model

A placement selects a set of distinct state vertices, one unit-cost dedicated
actuator per selected state. Any one selected actuator may fail; the state and
all state arcs remain. Structural controllability uses independent free
parameters. Write d for the maximum-matching deficiency of the state graph.

**Theorem.** It is NP-complete to decide whether a feasible instance has a
robust placement of size d+1. This holds for loopless bipartite digraphs with
exactly one source SCC, of size four. Every state is an allowed actuator
location, and no location is mandatory. In particular, all instances in this
reduction have positive deficiency and satisfy the absolute lower bound
OPT >= d+1.

Thus a polynomial exact algorithm for this model implies P=NP. More sharply,
the decision problem OPT <= d+1+k is para-NP-hard in k, already for k=0.
This does not contradict an FPT algorithm parameterized by d: d grows with
the Hamiltonian instance. It also supplies no constant inapproximability
factor: the NO/YES separation is only (d+2)/(d+1).

The proof below is direct. It does not rely on the existing sparse-hardness
or general-realization constructions, a fixed mandatory basis, or an
approximation-tight example.

## 1. Hamilton paths with specified ends

Start with a simple connected Hamiltonian Cycle instance G of minimum degree
at least two. This restriction preserves NP-completeness: disconnected graphs
or graphs with a vertex of degree less than two are immediate NO instances
and can be replaced by one fixed connected minimum-degree-two NO instance
(two triangles sharing one vertex).

Choose v in G. Replace v by two nonadjacent vertices s and t, each adjacent
to every former neighbor of v. Call the resulting graph H, with n vertices.
Then H has a Hamiltonian s-t path if and only if G has a Hamiltonian cycle:
split a cycle at v, or identify the two ends of the path back to v.
Here n >= 4, H is connected and has minimum degree at least two.

Create a resource for every ordinary edge of H, serving its two endpoints,
and two additional singleton resources a_s and a_t serving s and t. In the
pseudograph interpretation these singletons are loops at s and t. Let R be
the resource set and m=|R|. Initially the jobs are the n vertices of H.
All resources are vulnerable.

## 2. The state digraph

Use the following 2m+n+4 states:

- a resource state r_e and a center c_e for each e in R;
- one job state j_u for every u in V(H);
- four source states A, B, F_A, F_B.

The arcs are

    r_e -> c_e,
    c_e -> r_e and c_e -> j_u for every endpoint u served by e,
    A -> F_A -> B -> F_B -> A,
    F_B -> r_(a_s),       F_A -> r_(a_t),
    A -> c_e for every e in R.

No other arcs are present. The bipartition is

    X = {r_e} union {j_u} union {A,B},
    Y = {c_e} union {F_A,F_B}.

Every arc crosses the partition and there are no loops. The four source
states form an SCC with no incoming arc. Its outgoing star reaches every
center and hence every resource and job. No outside state reaches the source
cycle. It is therefore the unique source SCC. The construction is polynomial
and has O(|V(G)|+|E(G)|) states and arcs.

## 3. Its matching matroid

Let T be the transversal matroid on right-state vertices, using the left-state
copies as presentation rows, and M=T* its dual. With dedicated inputs,
a placement S satisfies the matching condition precisely when it spans M:
the unactuated states must be matchable by state arcs. Surviving every
selected-input failure is therefore robust spanning in M. Accessibility after
a failure is equivalent to selecting at least two source-SCC states.

The bipartition splits T into two independent summands. All Y columns are a
free summand: match c_e from r_e, F_A from A, and F_B from B. The extra star
arcs from A affect only this already-free summand. Thus the Y states are
loops of M and cannot improve matching rank.

On X, the rows are

    c_e: {r_e} union {j_u : u is an endpoint of e},
    F_B: {A,r_(a_s)},
    F_A: {B,r_(a_t)}.

They have full row rank m+2, matching each c_e to its private r_e and the two
filler rows to A,B. As |X|=m+n+2, rank(M)=n. Hence the constructed digraph has
matching deficiency d=n.

For clarity, remove A,B and the two filler rows temporarily. Call the
resulting dual matroid M_0. It is the assignment matroid on the original
resources together with one extra singleton resource b_u at every job u;
the column j_u corresponds to b_u. To verify this identification without
assuming a representation theorem, let Q be selected original resources and
J selected job columns. Duality and the private resource columns give

    r_(M_0)(Q union J) = |J| + nu(Q, V(H) minus J),

where nu is maximum matching between resources and their served jobs.
The same expression is the assignment rank after using the singleton
resources J first. An ordinary resource assigned to a job in J can always
be replaced there by its dedicated singleton, so this ordering loses no rank.

Adding a row {A,r_(a_s)} is a series extension on the transversal side:
contracting A recovers the former presentation. On the dual side it makes
A a parallel copy of r_(a_s), and similarly B a parallel copy of r_(a_t).
One can also check parallelism directly from the dual rank formula: deleting
both A and r_(a_s) loses precisely their filler row, while deleting either
alone does not lower full row rank. Their dual pair has rank one, and each
is a nonloop. Therefore M consists of M_0 with those two parallel copies,
plus the loops Y.

## 4. The absolute lower bound excludes the auxiliary choices

Any robust spanning set in a rank-n matroid has at least n+1 elements.
If a robust spanning set S has exactly n+1, deleting each element leaves
an n-element basis. Thus S itself is a spanning circuit. It contains no
matroid loop and, since n>=2, no parallel pair.

Suppose the constructed control graph has a robust placement S of size n+1.
The source quota requires two of A,B,F_A,F_B. The fillers are matroid loops,
so A and B must both be chosen. Their physical copies r_(a_s), r_(a_t)
cannot also be chosen. No center can be chosen either. Replacing A,B by
their corresponding original resources preserves size and all ranks.

We obtain a robust assignment of n+1 resources in the pseudograph H with
the two required endpoint loops a_s,a_t, and possibly the completion loops
b_u. This is why no cost-preserving domination lemma is needed here:
the absolute lower bound itself eliminates the extra singleton resources.

## 5. Bicyclic selections encode Hamilton paths

Take any robust assignment of n+1 pseudograph edges on the n jobs. Let a
connected component of its selected pseudograph have h vertices and q edges.
It has to cover all its h jobs after deletion of any selected resource in
that component, so q>=h+1. Isolated uncovered vertices are infeasible.
Summing over components gives n+1>=n+(number of components). Hence the
selected pseudograph is connected and has exactly two independent cycles.

Both prescribed endpoint loops are selected. Removing them leaves n-1 edges.
Connectivity therefore forces all these remaining edges to be ordinary
nonloop edges forming a spanning tree of H; no completion loop b_u can occur.
For every tree edge, failure of that resource splits the tree into two pieces.
Each piece must still match its vertex jobs. A tree piece has one fewer
ordinary edge than vertices, so it needs an endpoint loop. Consequently every
tree-edge cut separates s from t. Every tree edge lies on the unique s-t
path, so the tree is itself a Hamiltonian s-t path.

Conversely, a Hamiltonian s-t path together with a_s,a_t is a robust
assignment. If an endpoint loop fails, the remaining graph is a connected
unicyclic spanning graph. If a path edge fails, the two components are trees
with one endpoint loop each. In either case orient tree edges away from the
loop and assign every nonroot vertex its parent edge, assigning the root its
loop. This matches all jobs. Replace the endpoint resources by A,B to get a
robust control placement of size n+1; its source quota also holds.

This proves

    G Hamiltonian  <=>  H has a Hamiltonian s-t path
                   <=>  OPT_control = d+1,       d=n.

## 6. Feasibility, optional locations, and NP membership

All constructed instances are feasible. Every state has an incoming arc and
the only source SCC has four vertices. Selecting every state remains
structurally controllable after one input fails: use one state arc to match
the now-unactuated state and the surviving dedicated inputs for every other
state; accessibility is retained.

In fact no single location is mandatory. Every pair of distinct right-state columns
is matchable. Indeed each column either has at least two distinct predecessor
rows, or has a unique predecessor; no two degree-one columns share that
predecessor. This is immediate from the displayed construction (job columns
have at least two incident edge-resource centers). Thus T is simple.
For any prescribed state z, select all states other than z. After failure
of another selected input w, match the pair {z,w} using state arcs and all
other states by their dedicated inputs. The source cycle still has at least
two selected surviving actuators. This witnesses avoidability of every z.

A candidate placement is verifiable by polynomially many bipartite maximum
matchings, one per possible failed input, together with source accessibility.
The matching deficiency is polynomially computable. This establishes NP
membership and completes the proof.

## Scope and relationship to the existing results

- This strengthens the exact-complexity boundary: hardness already occurs at
  the theoretical minimum of one actuator above matching deficiency, with
  one four-state source SCC. It does not improve the existing constant
  approximation-hardness bounds.
- One source SCC does not mean that the whole digraph is strongly connected.
- This result excludes an FPT algorithm in excess k alone unless P=NP. It
  does not exclude FPT in matching deficiency, nor approximation schemes
  whose ratio does not distinguish d+1 from d+2.
- All graph states are available; there is no external mandatory-actuator
  constraint. A,B become necessary only for attaining the threshold.
- The standard Hamiltonian Cycle NP-completeness theorem is the only source
  hardness theorem used. The bicircular/assignment and matroid duality tools
  are classical. Priority of the particular control reduction has not been
  established; no first-result claim is made.

## Related literature

The graphic-matroid analogue of hardness at rank plus one is classical:
the graphic matroid of a connected v-vertex graph has rank v-1 and
admits a robust v-edge spanning set exactly when the graph has a Hamiltonian
cycle. Bentert, Fomin, Golovach and Morelle,
[Fault-Tolerant Matroid Bases, ESA 2025, Observation 3 and the preceding discussion](https://doi.org/10.4230/LIPIcs.ESA.2025.83),
recall this argument. It does not automatically yield the present control
result, since an arbitrary graphic matroid is not the dual state matching
matroid of a specified control instance, and source accessibility must also
be realized. The concrete cotransversal completion and single four-state
source gadget above establish that transfer directly.

[Adjiashvili, Bindewald and Michaels (2017)](https://arxiv.org/abs/1703.06074)
study the vulnerable-resource assignment model used as the intermediate
interpretation. Their general hardness results are not invoked as a black
box here. The source Hamiltonian Cycle theorem is classical
(Karp, *Reducibility Among Combinatorial Problems*, 1972).
Targeted literature searches do not establish publication priority for the
restricted control theorem; that question remains separate from its proof.

## Reproducible finite verification

Run `python verify_lower_bound_excess.py` with Python 3.10 or later. The
standard-library verifier passed:

- all 1,098 simple H graphs of orders two through five, inspecting 264,359
  anchored threshold supports including completion singleton resources;
- 38,846 actual state-plus-input threshold placements on 15 small control
  fixtures, with accessibility checked by direct graph reachability;
- five Hamilton-cycle-to-path-to-control fixtures, including K4, K3,3, the
  cube, the non-Hamiltonian Petersen graph, and two triangles sharing a
  vertex; all three positive Hamilton witnesses lift successfully;
- 4,333 omitted-state pairs, certifying individual location avoidability
  on those five targets, and 400 arbitrary placements compared against the
  explicit assignment representation.

The program enumerates candidate resources, matches actual input/state
bipartite graphs, and checks reachability separately from the proof's
matroid decomposition. Finite enumeration supplements the general proof.
