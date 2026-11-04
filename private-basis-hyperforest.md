# Private-basis reductions and the missing rank term

This extends `variable-basis-path-cover.md` to arbitrary option neighborhoods.
It gives an exact objective and obstructions to several Set Cover reductions;
it does not establish a stronger numerical hardness constant.

`assignment-hardness-transfer.md` proves that the same optimum formula
applies to every feasible unit-cost uniform robust assignment instance,
even without private resources: adding those resources preserves OPT by
a tight-Hall-set exchange.

Let U be n jobs. For each job j there is a private resource b_j. Every other
resource t_e can serve any job in a nonempty set S_e. Assume the S_e cover U;
resources with identical neighborhoods are still distinct elements.

## State realization and exact formula

Use states b_j, t_e, c_e, with bidirected edges c_e--t_e and c_e--b_j for
j in S_e. The matching deficiency is n. Its dual matching matroid is the
assignment matroid just described: complementing a basis exchanges the
matching of selected options to absent privates with the matching of
complementary targets to centers. This is the same basis-by-basis proof as
in the path-cover note, with S_e replacing an edge's two endpoints.
Centers are dual loops. Every incidence component has positive deficiency,
so robust spanning also supplies its source-SCC accessibility quota.
`private_realize` in `variable_basis.py` implements the construction.

Fix a set C of selected option resources, and write P for the selected
private jobs and F=U\P for omitted privates. For each j form the labeled
hyperedge N_C(j)={e in C : j in S_e} on vertex set C. Let r_C be the rank of
these n labeled hyperedges in the hypergraphic matroid. Singleton hyperedges
are matroid loops. Then

    OPT = n + min_{C covering U} (|C| - r_C).                 (1)

Here is a direct proof. Robust Hall is equivalent to

    |P intersect A| + |N_C(A)| >= |A|+1

for every nonempty A subset U. Necessity follows by deleting a selected
resource adjacent to A; sufficiency follows by Hall after any deletion.
It is equivalent to these two conditions:

1. C covers every job, including jobs whose private resource is selected.
2. For every nonempty A subset F, |N_C(A)| >= |A|+1.

For sufficiency, apply condition 2 to A intersect F when it is nonempty;
otherwise coverage gives the one extra neighbor. For necessity apply robust
Hall to a singleton job and to each subset of F. Condition 2 is precisely
the hyperforest independence criterion on the omitted job hyperedges.
Thus at most r_C privates can be omitted, and any maximum hyperforest
permits exactly that many omissions. The resulting cost is n+|C|-r_C.

The hypergraphic interpretation and its rank formula are classical. In the
notation here, [Frank--Király--Kriesell's rank formula, as stated in
Király, *Computing the minimum cut in hypergraphic matroids*, Theorem 1](https://egres.elte.hu/qp/egresqp-09-05.pdf)
gives the equivalent expression

    |C|-r_C = max_{partitions Pi of C} (|Pi| - cross_C(Pi)),  (2)

where cross_C(Pi) counts job hyperedges meeting at least two partition
blocks. The hyperforest criterion used above is also stated in
[Hörsch--Szigeti, *Packing of spanning mixed hyperarborescences with
flexible roots via matroid intersection*, Section 2](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v28i3p29/pdf/).
Hyperforests may equivalently be shrunk to ordinary forests, choosing two
vertices from each hyperedge. Consequently the excess over n is the minimum
number of components of such a forest on chosen option vertices that cover
all jobs; each connecting job may be used at most once.

Both solution directions for (1) are polynomial for a given C. From a
placement read C and its omitted privates. Conversely find a maximum
hyperforest greedily with an independence oracle: a job subset is independent
iff it remains matchable to C after deletion of each individual option.
These are polynomially many matching tests. Retain precisely the other
privates. This does not solve the outer minimization over C.

## Consequences for hardness constructions

The private basis is available but need not be mandatory. In fact if every
job has at least two options, the state transversal matroid is simple on
target columns: every pair has a two-row matching. Its dual therefore has
no two-element cocircuit, so no state is mandatory. This does not eliminate
the correction r_C in (1).

**Frequency at most two collapses.** If each job has at most two options,
the hyperedges on the full option set are ordinary edges or singletons.
When the incidence graph is connected their graphic rank is |C|-1.
Choosing all options in (1) gives OPT<=n+1; robustness in a rank-n matroid
gives the reverse inequality. For c incidence components OPT=n+c.
Thus deleting one forced private per covering requirement cannot preserve
an arbitrary Vertex Cover gap, regardless of the option set sizes.

**Connected regular incidences collapse too.** Suppose there are n options
and n jobs, and their connected bipartite incidence graph is k-regular,
k>=2. For every proper nonempty job subset A, regularity gives
|N(A)|>=|A|, and equality would disconnect A union N(A) from its complement.
Thus |N(A)|>=|A|+1. Selecting all n options and any one private satisfies
robust Hall, including A=U. Once again OPT=n+1. In particular connected
3-uniform, 3-regular Set Cover sources have no gap in this construction.

**Plain Set Cover can even order instances backwards.** On six jobs, two
disjoint triples have cover optimum two but robust optimum eight. In
contrast take the six edges of K4 as jobs and its four incident-edge triples
as options. Its cover optimum is three, whereas the frequency-two result
gives robust optimum seven. Selecting a larger cover can gain enough
hypergraphic rank to reduce the actuator cost. A gap proof must bound
|C|-r_C, not |C| alone.

For option neighborhoods of size at most three, all centers have state
degree at most four. Splitting high-degree job targets into parallel copies
as in `sparse_realize` gives maximum state degree four. The general
coloop-completion lemma in `hardness-obstructions.md` proves that these
parallel extensions preserve the optimum and gives a reverse solution map.
Odd edge subdivisions then provide the same polynomial girth padding as
in `hardness-degree-threshold.md`. This is available for a future hardness
reduction; it is not itself such a reduction.

## Exact verification

`python verify_private_basis_spares.py` compares (1) against direct resource
assignment ranks and all one-resource deletions for every covering multiset
of at most four option neighborhoods of size at most three on at most four
jobs: 2,200 instances. For every covering C it separately checks (2) by
enumerating all partitions and checks r_C by the strong Hall criterion on
every job subset. On sources through three jobs it also checks the actual
state-plus-input matching construction. It tests coloop completion on every
spanning resource support. Another 46 fixed-seed instances on five or six
jobs pass the same exact objective and partition checks.

The next unresolved source family has irregular incidences, option
neighborhoods of size at least three, and some jobs with at least three
options. A stronger lower bound would require a gap for the corrected
objective (1), or a different cotransversal construction.
