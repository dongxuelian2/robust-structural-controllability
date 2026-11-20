# Exact transfer of the Max-E3-Lin-2 gap

**Theorem.** Approximating unit-cost 1-RMSCP-d is NP-hard within
$843/838-\epsilon$ on general digraphs and within
$1474/1469-\epsilon$ on strongly connected digraphs, for every fixed
$\epsilon>0$. Both bounds hold for loopless bipartite state graphs with
underlying degrees in $\{2,3\}$, maximum indegree and outdegree three,
no mandatory actuator location, and girth at least $N^{1-\delta}$ for
every fixed $0<\delta<1$. The strong targets are bidirected.

The construction uses the shared backup from
[the Vertex Cover reduction](shared-backup-hardness.md).
It works directly on amplified parity equations. Its cost in excess
of the baseline is exactly the minimum number of unsatisfied equations.

## Literal choices and private guards

Start with $q$ equations in three variables over $\mathbb F_2$. For
each variable use a degree-two/three amplifier as in Definition 7 of
Chlebík–Chlebíková. A contact has degree two and is assigned to one
occurrence of that variable. All occurrences have distinct contacts.
A checker has degree three. In each amplifier, the number of edges
across any cut is at least the number of contacts on the side with
fewer contacts. Write $M$ for the total number of amplifier vertices.
There are $3q$ contacts.

Each amplifier vertex has a pair of resources, one for each Boolean
value. Make a job for the pair and give all $M$ such jobs a single
shared backup. A guard job served only by that backup and its guard-only
partner forces both resources. On backup failure each pair job needs
at least one of its two labels.

For each amplifier edge, use the two disjunctive clauses expressing
equality of the endpoint labels. For each parity equation, use the four
three-literal clauses excluding its four falsifying assignments. Each
clause has its own job, served by the labels making any of its literals
true and by a private resource. In the general variant force the
private resources with paired quota-two groups. Use the guard-only
partner to complete the last pair if needed. If the private count is
even, put the shared backup and its partner in their own group. In the
strong variant use a separate two-resource guard for each private
resource, as in the Vertex Cover construction.

When the private resource for a clause fails, its job needs at least
one satisfying selected label. In the strong variant the same is true
when either of its guard resources fails. Conversely, if each pair and
each clause has a label, all jobs survive one failure. Normally each
pair job uses its own label, and clause jobs use their private resources.
A failed label affects at most one pair job. A label diverted to repair
a failed private resource also affects at most one pair job. In both
cases the shared backup takes that job if needed. If one of the two
shared guard resources fails, the other takes the guard and the
ordinary label assignments suffice. Hence feasibility is exactly the
pair and clause condition on the chosen labels, plus the forced
resources.

The amplifiers have $(3M-3q)/2$ edges. The number of clause jobs is
$p=2(3M-3q)/2+4q=3M+q$. For a feasible selection with $D$ doubly
selected pairs, the two costs are $M+D+p+2$ and $M+D+2p+2$.
An assignment violating $u$ original equations gives $D=u$. Copy
each variable's value throughout its amplifier, and add the other
label at one contact of each violated equation. This satisfies its
four disjunctions, and preserves all the equality clauses.

## Rounding one violation per doubly selected pair

Given a feasible selection, assign all singly selected vertices their
fixed labels and assign doubly selected vertices arbitrary labels.
Let $K$ count disagreeing amplifier edges and let $U$ count unsatisfied
parity equations on the contact labels. Every disagreeing edge has at
least one doubly selected endpoint. Every unsatisfied equation has at
least one doubly selected contact.

Repeatedly flip a doubly selected vertex if this strictly decreases
$K+U$. This terminates in polynomial time, since the potential is an
integer bounded by the number of edges and equations. At termination:

- A doubly selected checker meets at most one disagreeing edge. If
  two or three incident edges disagreed, flipping would lower $K$.
- A doubly selected contact meets at most one disagreeing edge.
  Flipping when both edges disagree lowers $K$ by two and changes
  $U$ by one in absolute value.
- A doubly selected contact of an unsatisfied equation meets no
  disagreeing edge. If it met exactly one, flipping would keep $K$
  fixed and fix that equation.

Charge every disagreeing edge to one doubly selected endpoint.
Charge every unsatisfied equation to one doubly selected contact.
These charges are all distinct by the three properties and the
one-equation incidence of each contact. Thus $K+U\le D$.

Now give each original variable the majority label among its contacts.
By the amplifier cut inequality, at most $K$ contacts change in total.
Each change affects one equation, so the resulting assignment violates
at most $K+U\le D$ equations. With $u^*$ the source optimum number of
violations, the exact costs are therefore

$$
 \operatorname{OPT}_{\mathrm{gen}}=4M+q+2+u^*,\qquad
 \operatorname{OPT}_{\mathrm{sc}}=7M+2q+2+u^* .
$$

## Realization and the constants

The full set of original resources is feasible. The assignment and
quota realizations from the shared-backup reduction apply unchanged.
In reverse, remove all added private and branching resources by the
dominated-resource exchange, then round the retained original labels.
The exchange preserves selected originals and the pair quotas. The
matching rank exceeds the number of forced resources by $M-1$ in
both variants. Every individual label is avoidable in the abstract
instance by selecting the other label and both labels at all other
vertices.

Use the balanced hardness source in Chlebík–Chlebíková, Theorem 2.
Each variable has $3k$ occurrences for a sufficiently large fixed $k$.
Their amplifier construction uses at most $(207/10)k$ vertices for
$3k$ contacts, so $M\le(207/10)q$. The fixed amplifier can be
hard-coded for each chosen source gap. The source gap has $u^*\le\eta q$
in the good case and $u^*\ge(1/2-\eta)q$ in the bad case.
Disjoint replication with a single shared backup absorbs the additive
two. Taking $\eta$ arbitrarily small gives the limiting ratios

$$
 1+\frac{1/2}{4(207/10)+1}=\frac{843}{838},\qquad
 1+\frac{1/2}{7(207/10)+2}=\frac{1474}{1469}.
$$

The source gap and amplifier constants are from
[Chlebík–Chlebíková, Theorem 2 and Section 3](https://pure.port.ac.uk/ws/portalfiles/portal/1887750/3DM_JOURNAL_revision_old.pdf).
`verify_direct_lin2.py` checks both cost maps with exact source
matchings and checks actual actuator failures and recovery on two
padded control graphs. The $3/2-\epsilon$ target remains open.
