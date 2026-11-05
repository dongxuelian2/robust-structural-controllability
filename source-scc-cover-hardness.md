# Using source SCCs to halve the covering baseline

The unit-cost, single-actuator-failure model on general directed graphs
admits a Set Cover reduction with normalized cost |U|+tau, rather than
2|U|+tau. The exact identity, including a two-copy parity construction, is

    OPT(G) = 2*(|U| + tau).

All states remain selectable. For a frequency-two covering source, G is
loopless and bipartite, with every underlying undirected degree in {2,3}
and both maximum indegree and maximum outdegree at most three. It is weakly
connected but has several source SCCs. This gives **250/249-epsilon**
NP-hardness of approximation on this class. The previously proved
**2238/2233-epsilon** bound still applies to the separate, strongly
connected class in `hardness-degree-threshold.md`.

The change is to force a private basis through accessibility, pairing its
members in source two-cycles. They need not be mandatory for matching.
Consequently the forced-matching-basis obstruction in
`hardness-obstructions.md` does not apply.

## Construction and both solution maps

Let nonempty sets S_1,...,S_q cover a nonempty universe U of size n. Take
two disjoint copies indexed by c=0,1. Introduce:

- a private-basis state b_(c,j) for each job j in each copy;
- an option state a_(c,i) for each set;
- a physical job state z_(c,j).

Pair b_(0,j),b_(1,j) in a source directed two-cycle. Put all option states
in a single directed cycle, alternating copies in the order
a_(0,1),a_(1,1),...,a_(0,q),a_(1,q). These are the source SCCs. Job states
have no outgoing arcs.

For a source-cycle state v, write p(v) for its predecessor. Add
p(b_(c,j))->z_(c,j). The option p(a_(c,i)) reaches precisely the jobs
z_(c,j) with j in S_i through its own binary tree. Use one outgoing arc
from the source cycle into that tree. An internal tree node is a pair
x,y with x->y, y->x; its incoming parent arc ends at x and its at most
two child arcs start at y. A singleton option can connect directly to its
job. All trees are disjoint except for shared job leaves.

Every private-basis pair is a source SCC of size two, so all 2n states

    B = {b_(c,j)}

are mandatory for robustness. The option cycle requires two selected
actuators. The remaining SCCs are tree pairs and job sinks, all with
incoming arcs; they impose no additional source quota.

**Matching.** Match every source-cycle column from its predecessor row
and every internal pair across its two-cycle. This saturates all nonjob
columns. There are exactly 2n job states and their rows are empty, so
the matching deficiency is exactly 2n. Actuating B allows its predecessor
rows to match all jobs privately; B is a basis of the dual matching
matroid M.

Fundamental exchanges from B have a direct description. An option column
can replace precisely the b_(c,j) for j in its set. A tree state x can
replace precisely the basis states for leaves below it, while y is an
M-loop. A physical job z_(c,j) replaces its corresponding basis state.
For example, selecting an option frees its predecessor row; alternating
through the tree's matched pairs frees the row needed for any one of its
job leaves. No path reaches another option's tree or the other copy.

Equivalently, permute the source-cycle rows so their distinguished cycle
columns become diagonal. An internal row x has only column y; that column
is a free summand of the transversal matroid and can be removed with its
row. The remaining internal row has {x} and its child columns. This is
exactly the disjoint resource-tree presentation from
`assignment-hardness-transfer.md`, together with the private basis B.
It proves both the positive and negative exchange assertions, including
the absence of other fundamental neighbors.

Since B is selected, failures outside B leave a spanning basis. Robust
matching is therefore equivalent to the extra selected states covering
all 2n basis elements by their fundamental neighborhoods.

Given two source covers C_0,C_1, select B and their option states. They
cover all basis elements, so matching survives every failure. Each cover
is nonempty, hence the option cycle has at least two actuators. After any
failure every source SCC still has an actuator. This is a feasible
placement of size 2n+|C_0|+|C_1|.

Conversely, take any robust placement S. It contains B. Map every selected
option or functional tree state to its owning set in its copy. Map a
selected physical job to any source set containing that job, and discard
selected M-loops. Its fundamental-neighborhood coverage shows that the
resulting C_0,C_1 each cover U, with

    |C_0|+|C_1| <= |S|-2n.

Taking the smaller of these two covers gives a source solution of size at
most |S|/2-n. Both maps are polynomial, proving the exact OPT identity
and its approximation transfer. Tree states and physical job inputs
cannot lower the optimum.

**Degree and parity.** Give every b_(c,j), a_(c,i), z_(c,j), and tree x
the color c; give its tree partner y color 1-c. Every arc crosses colors,
including both source-cycle types. Thus the underlying graph is
bipartite. Source basis vertices have undirected degree two; option
vertices have degree at most three because only one tree arc leaves a
cycle vertex. Tree x has degree two and y degree three. Each job has one
private incoming arc and one leaf arc per source set containing it. At
frequency two it has degree three. In general its degree is frequency
plus one, and every outdegree is still at most three. There are no loops
or leaves. The option cycle and the covering incidences make the graph
weakly connected.

The state transversal matroid is simple: each column has at least one
incoming neighbor; columns with just one have distinct such neighbors;
all other columns have at least two. Thus every pair of columns matches.
No state is mandatory for robust matching alone. B is mandatory through
accessibility, not through a two-element cocircuit.

## The numerical gap

Use Vertex Cover on a 4-regular graph H with v vertices and m=2v edges.
Its edge universe and incident-edge vertex sets have frequency two.
The construction gives

    OPT(G(H)) = 2*(m + tau(H)).

For maximum source degree Delta, m<=Delta*tau(H), and an alpha-approximate
placement gives a source cover of size at most

    alpha*tau(H) + (alpha-1)*m
       <= [1+(Delta+1)*(alpha-1)]*tau(H).

The stronger absolute source gap is the one already audited in
`hardness-degree-threshold.md`: [Chlebik--Chlebikova, Theorem 17 and
Section 3.2](https://pure.port.ac.uk/ws/portalfiles/portal/1887750/3DM_JOURNAL_revision_old.pdf)
give, for arbitrarily small fixed eta>0, NP-hardness of distinguishing

    tau(H) <= v*(2*mu+8+2*eta)/(4*mu+12)

from

    tau(H) >= v*(2*mu+9-2*eta)/(4*mu+12),

where mu<=217/10+eta. The common factor two in the target cancels.
Its NO/YES ratio is

    1+(1-4*eta)/(10*mu+32+2*eta)
       >= 1+(1-4*eta)/(249+12*eta).

It tends to 250/249 as eta tends to zero. Choosing eta sufficiently small
proves 250/249-epsilon hardness for every fixed epsilon>0 for which the
factor exceeds one, unless P=NP. Only the source theorem, not a new PCP
claim, enters this calculation.

[Pocas--Pequito--Ramos, Automatica 190 (2026), 113074](https://doi.org/10.1016/j.automatica.2026.113074)
use OPT=2|U|+tau. Their construction and the uniform robust-assignment
constructions cited in `hardness-degree-threshold.md` force the baseline
through matching. Here source SCCs instead force each private basis
element once. This targeted comparison establishes the change from those
constructions; priority beyond the sources reviewed is unconfirmed.

This does not give 3/2 hardness. Even this smaller baseline still dilutes
the source gap, and these graphs are not strongly connected. Merging
their source SCCs would remove the constraints that force B, so the
strongly connected realization theorem cannot simply be applied to
retain this objective.

## Exact certificates

`scc_cover.py` builds both the binary-tree construction and a direct-port
version, and implements both cover maps. `verify_scc_cover.py` checks
the actual state-plus-input matching after failures of every forced
actuator. It compares all placements containing B when there are at most
12 other states, and fixed-seed placements otherwise, against the full
predicted coverage and source-quota conditions. It additionally verifies
the exact source SCCs, deficiency, simplicity test, the OPT identity on
the exhaustive instances, and the recovered cover maps. The K5, four-cycle
and K4 Vertex Cover fixtures check the loopless, bipartite subcubic
construction, both directed degree bounds, and explicit optimum witnesses.
