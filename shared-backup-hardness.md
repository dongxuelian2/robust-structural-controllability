# Hardness from sharing a backup across a matching

**Theorem.** For every fixed $\epsilon>0$, approximating unit-cost
1-RMSCP-d within $913/908-\epsilon$ is NP-hard on general digraphs.
The strongly connected restriction is NP-hard within $1594/1589-\epsilon$.
Both bounds hold on loopless bipartite targets with underlying degrees
in $\{2,3\}$, maximum indegree and outdegree three, and no mandatory
actuator location. The strongly connected targets are bidirected.
For any fixed $0<\delta<1$ both reductions permit undirected girth at
least $N^{1-\delta}$ on $N$ states.

The reductions take a Vertex Cover instance with a designated matching
and vertex-disjoint triangles. A single backup suffices for all matching
edges because one unavailable vertex affects at most one matching edge.
The resulting cost depends on the edges left after removing both structures.

## A shared backup and triangle quotas

Let $H=(V,E)$ be a simple graph with vertex-disjoint triangles $\mathcal T$
and a nonempty matching $P$ containing no edge internal to these triangles.
Let $R$ be the edges outside $P$ and the triangles. Use one job per edge of
$P\cup R$. Vertex resource $v$ serves any of its incident jobs. Each job
$e\in R$ has its own private resource $b_e$. All jobs of $P$ share one
backup $b$. A guard job is served only by $b$ and a guard-only resource $g$,
so both are forced by robustness.

Pair the $b_e$ into quota-two groups, using $g$ to complete the last
pair if needed. When their count is even, put $b,g$ in their own group.
For each triangle impose quota two on its three vertex resources.
Assume every vertex is incident to a retained edge, as holds for the
hard family below.

Any feasible selection contains the forced resources. On failure of
$b_e$ for $e\in R$, its job needs a selected endpoint. On failure of $b$,
every $P$ job needs a selected endpoint. The triangle quotas ensure that
their internal edges are covered. Thus the selected vertices form a
vertex cover.

Conversely, select any vertex cover and all forced resources. With all
selected vertices available, assign distinct cover vertices to the $P$
jobs by the matching property. If a vertex fails, at most one such job
loses its resource and $b$ fills in. If a private $b_e$ fails, give its
job a selected endpoint and use $b$ for the at most one $P$ job displaced.
In both cases $g$ takes the guard job. On loss of $b$ or $g$, the other
takes the guard and all edge jobs retain their ordinary assignments.
Thus with $\tau(H)$ the minimum cover size,

$$
 \operatorname{OPT}_{\mathrm{gen}}=|R|+2+\tau(H).
$$

## Enforcing the guards and triangles by matching alone

For the strongly connected variant, use the same shared backup and
guard. For each $e\in R$ give $b_e$ its own guard job, served only by
$b_e$ and a second private resource that serves just this guard.
Both are forced. At an ordinary assignment they can serve the edge
job and its guard. On loss of either one, the survivor must serve the
guard and the edge job requires a selected endpoint.

Replace each triangle quota by a job served only by its three vertex
resources. Robustness requires at least two of them to be selected.
All other soundness arguments are unchanged, so every robust selection
costs at least $2|R|+2+\tau(H)$.

For completeness, it suffices to lift a cover with the following local
property. The packing edges incident to a triangle form either an empty
set or two edges from distinct triangle vertices to two outside contacts,
of which at least one is selected. The contact pairs of different
triangles are disjoint and use no triangle vertices. The hard family uses
exactly two such edges at every triangle.

In the ordinary assignment use a selected contact when possible. At
most one packed job per triangle then needs a triangle vertex, leaving
another selected vertex for the triangle job. If a selected resource
fails or a selected vertex must be diverted to a job in $R$, at most one
of these ordinary assignments is lost. If it is a $P$ job, use $b$.
If it is a triangle job and another selected vertex of its triangle
was assigned to a $P$ job, move that vertex to the triangle and use
$b$ at the displaced job. Otherwise use a spare selected vertex of
the triangle. Jobs in $R$ and the guards are handled as above.
Thus any such cover lifts at cost $2|R|+2+|C|$.

## The hard source and two packings

Use the degree-five construction of Chlebík–Chlebíková, Sections 2,
3.3 and 4, in the proof of Theorem 17. They reduce a balanced
Max-E3-Lin-2 instance with $q$ equations to a graph with cover optimum

$$
 \tau(H)=(\mu+5-\beta)q.
$$

Here $\beta$ is the optimum satisfied fraction and $\mu$ is the number
of amplifier vertices per $k$ equations at each variable. Their
amplifiers allow $\mu\le207/10$ for sufficiently large fixed $k$.
Their standard independent sets give minimum covers that select
exactly one vertex of each contact pair. Their gap for $\beta$ is
$1-\eta$ versus $1/2+\eta$ for arbitrarily small fixed $\eta>0$.
[Chlebík–Chlebíková, Theorems 2, 16 and 17](https://pure.port.ac.uk/ws/portalfiles/portal/1887750/3DM_JOURNAL_revision_old.pdf).

For the degree-five construction each amplifier vertex is doubled.
The two copies are adjacent except at the $k$ contacts of the third
position. Each amplifier edge creates two cross edges. The $3k$ contacts
have degree two in the amplifier and all checkers have degree three.
An amplifier on $\mu k$ vertices therefore contributes
$2\mu k$ vertices and
$(\mu k-k)+(3\mu k-3k)=(4\mu-4)k$ edges.
The equation gadget $G_j[4]$ contributes six additional vertices and
seventeen edges. Hence $H$ has $(2\mu+6)q$ vertices and $(4\mu+13)q$ edges.

In $G_0[4]$ write $a_0,a_1$ for the accented $x$ vertices and $s_{ijk}$
for the four special vertices. The triangles are

$$
 \{a_0,s_{101},s_{110}\},\qquad \{a_1,s_{000},s_{011}\}.
$$

The following four edges occur in the published gadget:

$$
 y_0s_{011},\qquad y_1s_{000},\qquad
 z_0s_{101},\qquad z_1s_{110}.
$$

These adjacencies and the two triangles can also be read from the edge
list in the proof of Theorem 19(B). For $G_1[4]$ complement the bit indices.
Contacts are distinct across equations, so all equation triangles
are vertex-disjoint.

- **General.** Use the four displayed edges and $x_1a_0,x_0a_1$ in $P$
  at each equation. Match all remaining vertices by pairing the two
  copies of each checker. This is a perfect matching disjoint from
  the triangle edges. Thus $|P|=(\mu+3)q$ and $|R|=(3\mu+4)q$.
- **Strongly connected.** Use the four displayed edges and the edge
  $x_0x_1$ from the consistency gadget. Again pair the two copies of each
  checker. Now $|P|=(\mu+2)q$ and $|R|=(3\mu+5)q$. The two contacts of
  each triangle are the two values of the same variable occurrence.
  A standard minimum cover selects one of them, as required for
  completeness. Thus the strong construction attains its lower bound.

The two exact optima on these graphs are

$$
 \operatorname{OPT}_{\mathrm{gen}}=(4\mu+9-\beta)q+2,\qquad
 \operatorname{OPT}_{\mathrm{sc}}=(7\mu+15-\beta)q+2.
$$

To absorb the constant, take disjoint copies of the source graph and
keep a single shared backup for their union of matchings. All source
counts scale except the shared backup and its guard partner.
As $\eta$ tends to zero and $q$ grows, the two NO/YES ratios are bounded
below in the limit by

$$
 1+\frac{1/2}{4(207/10)+8}=\frac{913}{908},\qquad
 1+\frac{1/2}{7(207/10)+14}=\frac{1594}{1589}.
$$

## Control realization and checks

Use [uniform assignment transfer](assignment-hardness-transfer.md) for
the strong variant. For the general variant apply the same cotransversal
presentation followed by [exact quota realization](general-sparse-equivalence.md).
The branching resources and private job resources are dominated by
original resources. The reverse exchange retains selected originals,
so neither the forced pairs nor the triangle quotas can be bypassed.
Parallel and loop additions, leaf completion and girth padding preserve
the cost bounds and keep all locations selectable.

The forced resources alone lack matching rank in both constructions.
In the general variant the job rank is $|R|+|P|+1$, whereas the forced
set has size $|R|+2$. In the strong variant the corresponding counts
are $2|R|+|P|+|\mathcal T|+1$ and $2|R|+2$. In the hard family a full
selection minus any one vertex resource is feasible. Thus there is no
fixed mandatory basis, and cover choices supply the missing rank.

`verify_shared_backup_cover.py` enumerates the supports of both
equation gadgets for both constructions. It also checks direct
state-plus-input matching under all selected actuator failures in
two padded control fixtures, and recovers covers after selecting
additional auxiliary states. The two constants improve the earlier
[general](source-scc-cover-hardness.md) and
[strongly connected](hardness-degree-threshold.md) reductions.
The $3/2-\epsilon$ target remains open.
