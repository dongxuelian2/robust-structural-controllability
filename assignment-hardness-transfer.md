# Arbitrary uniform robust assignments have subcubic control realizations

Let U be n jobs and E a set of distinct unit-cost resources. Resource e can
serve any one job in its nonempty neighborhood S_e. Every resource is
vulnerable. A selected set is robust if, after any one selected resource
fails, the jobs can still be matched injectively to surviving resources.
Assume the full original resource set E is robust. This is checkable by
polynomially many matching tests.

**Theorem.** Every such instance has a polynomial-size, loopless, bidirected,
bipartite, strongly connected control realization of maximum degree three
with exactly the same optimum. Both solution maps are polynomial and do
not increase cost. There is no additive actuator baseline. For each fixed
0<delta<1 the realization can also have undirected girth at least
N^(1-delta), where N is its number of states. Thus every approximation gap
for unit-cost uniform robust assignment transfers with its full ratio,
without any restriction on the original resource degrees.

The proof uses the sparse cotransversal realization in
`cotransversal-hardness-transfer.md` and the following exchange lemma. It
does not claim a new numerical lower bound or literature priority for the
exchange lemma. The source model is the uniform cardinality case of
[Adjiashvili--Bindewald--Michaels, *Robust Assignments with Vulnerable
Nodes*](https://arxiv.org/abs/1703.06074); extending the source model to
arbitrary costs or specified vulnerable subsets is not part of the theorem.

## Dominated added resources do not improve the optimum

**Lemma.** Add any number of resources whose neighborhoods are each
contained in the neighborhood of some original resource. Every robust
selected set of the enlarged instance can be converted to a robust set of
original resources of no larger cardinality, in polynomial time.

**Proof.** Robust Hall says

    h(A) := |N_S(A)|-|A| >= 1       for all nonempty A subset U.

Consider one selected added resource a, with neighborhood Q contained in
the neighborhood of original resource e. If e is unselected, replace a by
e. Suppose e is selected. Call A tight if h(A)=1, and consider the family
of tight A meeting Q. If there are none, deleting a is safe.

Otherwise this family is closed under intersection, with each intersection
still nonempty and meeting Q. Indeed two members A,B have both a and e
as common selected neighbors. If A,B were disjoint, their union would
satisfy h(A union B)<=0. If they intersect, submodularity and robust Hall
give h(A intersect B)=h(A union B)=1, with equality also in neighborhood
submodularity:

    N_S(A) intersect N_S(B) = N_S(A intersect B).

In particular a is a neighbor of the intersection, so it meets Q. Taking
the intersection of the finite family therefore gives a nonempty tight
set K contained in every tight set meeting Q.

The selected neighbors of K number |K|+1 and include a. Thus at most |K|
original neighbors of K are selected. Feasibility of the original instance
gives |N_E(K)|>=|K|+1, so an unselected original neighbor w of K exists.
Replace a by w. Every tight set losing a gains w, and every other job set
loses at most one unit of Hall slack. Robustness is preserved. Repeat for
all selected added resources. To implement the exchange, try deleting a,
then try each unselected original replacement, using matching to test
robustness. The proof guarantees success; computing K is unnecessary.

Identical resource neighborhoods are allowed in this lemma. They need not
be parallel elements of the assignment matroid: two copies of a resource
serving two jobs can both be used in one matching. This Hall exchange is
distinct from the matroid parallel-extension lemma.

## Binary branching and the exact control reduction

First add a private resource b_j serving each job j. For each original
resource e, build a directed binary tree with root t_e and one leaf b_j for
every j in S_e. The trees have disjoint internal vertices and share only
job leaves. Direct their arcs away from their roots. Jobs are the n sinks;
every vertex, including every branch vertex and sink, is selectable. With
k=|S_e| leaves, use max(1,k-1) internal vertices. A chain of binary forks
suffices; balanced trees are unnecessary.

Let M be the strict gammoid linking selected vertices to these sinks. Its
rank is n. On the original resource roots and the private job vertices,
its restriction is exactly the augmented assignment matroid: every
assignment uses private, mutually disjoint tree routes, and every linking
assigns its distinct terminal jobs to its distinct starting resources.
In particular every original robust assignment remains robust in M.

For the reverse direction, relax vertex-disjointness inside the trees.
This gives an assignment instance R on the same ground set: a vertex can
serve every descendant job. Every M-independent set is R-independent.
Both matroids have rank n, so an M-robust set is R-robust. Each branch
resource in R is dominated by its original root resource; each private
job resource is dominated by any original resource serving that job.
The lemma therefore converts an M-robust set to an original robust
assignment of no larger cardinality. Hence rho(M)=OPT_assignment.

For each internal vertex v use the transversal row

    {v} union out-neighbors(v).

These rows present M*, by the standard strict-gammoid/transversal duality.
Each has at most three elements. The rows have full rank by their distinct
identity entries v, and no element is a loop: internal vertices have their
identity entry, and every job has an original resource serving it. If I is
the number of internal vertices, there are n+I elements and I rows. Apply
the sparse cotransversal theorem to obtain the claimed control graph.
Its additional targets are parallel copies and its centers are dual loops;
coloop completion removes these copies before the Hall exchange removes
the branch and private resources.

Writing q=sum_e |S_e|, we have I<=q and the unpadded control graph has

    N <= 10I-n <= 10q-n.

Its matching deficiency is n. The graph is strongly connected, so its
accessibility quota is automatic for robust sets, which already have at
least n+1>=2 elements. The same odd-subdivision padding as in the sparse
cotransversal theorem supplies the girth bound and preserves both maps.
If every original resource can be avoided by some original robust set,
no control actuator is mandatory: the original solutions avoid every new
branch resource, and they lift while avoiding any prescribed copy of an
avoidable original resource.

For an unconditional version without mandatory locations, compose the
construction with `complete_leaves` as proved in
`row-splitting-equivalence.md`. This gives all degrees in {2,3}, preserves
OPT and the same reverse map, and permits the same girth padding.

The reduction is implemented in `assignment_transfer.py`. A source set
lifts to its retained original-root indices. To reverse a placement, call
`assignment_from_placement`, including after `subdivide_ports` padding.
All actuators remain distinct, of unit cost, and vulnerable to failure.

## The hyperforest formula applies without an available private basis

The lemma shows that adding one private resource per job preserves the
optimum of every feasible unit-cost uniform robust assignment instance.
Consequently `private-basis-hyperforest.md` gives, for the original instance,

    OPT = n + min_{C subset E covering U} (|C|-r_C),

where r_C is the hypergraphic rank of the job hyperedges N_C(j) on the
chosen resource vertices C. This formula does not require private resources
to be present in the source. A stronger assignment hardness construction
must control this rank correction even if it omits private resources.
The outer minimization remains; the identity is not a polynomial algorithm.

## Exact checks

`python verify_assignment_transfer.py` exhausts all feasible multisets of
at most six nonempty resource neighborhoods on up to three jobs: 979
sources. It checks every robust set using one added dominated resource,
128,128 exchanges in total, against direct assignment ranks. It separately
checks the branching rank inequalities, original-resource restrictions,
normalization, actual state-plus-input matching after each failure, and
reverse maps from arbitrary selectable control vertices. These checks
exhaust feasible sources on up to three jobs and four resources and include
larger, irregular and padded examples. One disconnected four-job source
has optimum six, strictly above the rank-plus-one lower bound.
