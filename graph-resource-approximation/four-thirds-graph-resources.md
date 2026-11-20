# A 4/3 approximation for graph resources with mandatory singletons

## Scope

This is a strict improvement over 3/2 on the simple graph-resource subclass with mandatory or forbidden singletons, including the strongly connected path family that makes arbitrary fixed-basis repair asymptotically 3/2, and the full control instances in the absolute-lower-bound Hamiltonian hardness reduction. It is not a 4/3 algorithm for arbitrary cotransversal matroids or arbitrary state graphs. No priority claim is made.

Let H=(V,E) be a simple undirected graph. An ordinary resource e serves its two endpoints. Each vertex may also have any number of distinct unit-cost singleton resources serving that vertex only; depict these as loops. A subset F of singleton resources may be mandatory; other singleton resources may be forbidden and removed in advance. Select a minimum-cardinality resource set whose incidence bipartite graph can match every vertex after deletion of any selected resource. All resources are vulnerable. The theorem allows arbitrary numbers of singleton copies; it does not allow parallel ordinary edges.

**Theorem.** This problem has a polynomial-time 4/3 approximation, including mandatory singleton resources. A fixed nonnegative additive charge (for selected matroid-loop states) preserves the ratio.

## 1. Exact feasibility characterization

For a selected pseudograph S, robust assignment is equivalent to the following two conditions:

1. Every vertex is incident to at least two distinct selected resources (a singleton counts once here).
2. No connected component of S is an ordinary simple cycle.

Equivalently, using standard graph degree (loops count twice), every vertex has degree at least two and no connected component is unicyclic.

Proof of the latter formulation. Necessity of graph degree at least two follows since a vertex of degree zero or one cannot survive its sole incident resource's failure. In any connected component with h vertices and q selected resources, robustness requires q>=h+1, since deleting any resource must leave enough resources to cover its h jobs. Thus no component can be unicyclic.

Conversely suppose every degree is at least two and every component has at least two cycles. Delete one resource. If it is not a bridge (including when it is a loop), its component stays connected and has at least as many edges as vertices, so contains a spanning unicyclic subgraph. If it is a bridge, neither resulting component can be a tree: such a tree would have at least two leaves, but only the endpoint of the removed bridge could have lost degree; a singleton resulting component would have had original degree one. Thus both pieces have cycles and spanning unicyclic subgraphs. Assign each unicycle's cycle resources bijectively to its cycle vertices, then assign every remaining vertex its parent tree edge. This covers all jobs.

For the equivalence with the first formulation, a finite connected unicyclic graph with minimum graph degree two is precisely a cycle, possibly a loop or a double edge. H is simple, so a double edge is unavailable. A single-loop component violates the two-distinct-resources condition. A component containing a loop and an ordinary edge has a vertex of degree at least three if every vertex has two incident resources, hence cannot be unicyclic. These observations establish the stated equivalence.

## 2. Preprocessing and canonical singleton resources

Work separately in every connected component of H.

- At an isolated vertex, feasibility requires at least two available singleton resources. Select all mandatory resources and, if necessary, add arbitrary available ones until two are selected. This is exact.
- At an ordinary-degree-one vertex, at least one singleton resource is necessary. If none is available, the instance is infeasible. If no singleton there is mandatory, designate one available singleton as retained. Unit costs make its identity immaterial: replace any chosen singleton at that vertex by the designated one when needed.
- If an ambient component is an ordinary simple cycle and contains no mandatory singleton, its optimum is exactly |V(component)|+1: the general resource-count lower bound applies, and all its cycle edges plus any available singleton attain it. If no singleton is available, it is infeasible. Solve this component exactly and remove it.

On every remaining component let K be all mandatory singleton resources, together with the designated singleton at each leaf not already having a mandatory singleton. There exists an optimum solution containing exactly K and no other singleton resources.

Here is the full elimination argument. Start with an optimum solution and ensure it contains each designated leaf singleton by an equal-cost swap at that vertex. Consider a selected singleton ell outside K, located at v. Remove ell temporarily.

If the result is still robust, keep the deletion. Otherwise its only potentially damaged component is ell's former component. This component was connected and remains connected after loop deletion. If its graph minimum degree drops below two, v has no remaining loop and at most one selected ordinary incident edge. It cannot be an ordinary leaf, because a retained leaf loop would still be present, and it is not isolated. Therefore v has an unused ordinary incident edge in H. Add that edge in place of ell. It restores the degree at v, does not lower any other degree, and either increases the cycle count of the damaged component or joins it to another previously robust component. The result is robust.

If all degrees remain at least two but robustness fails, the damaged component is unicyclic and hence a cycle. It is either an ordinary simple cycle or a lone singleton loop. In the latter case, v is not isolated, so an unused ordinary edge touches it. In the former case, an unused ordinary edge touches its vertex set unless the entire ambient component H is exactly that cycle. That exceptional ambient component either was removed in preprocessing or has a mandatory singleton; the latter would remain selected on the cycle and contradict that the damaged component is an ordinary cycle. Thus in every case there is an unused ordinary edge touching the damaged component. Adding it increases the cycle count if it is internal, or joins the damaged unicycle to another robust component if it is external. It restores robustness.

Each operation deletes one nonretained singleton, adds at most one ordinary edge, and introduces no singleton. Iterate. The resulting robust solution contains exactly K and costs no more. In fact if the starting solution was optimal it has the same optimum cost.

## 3. Polynomial relaxation and augmentation

Write k_v for the number of retained singleton resources at v. Set

- b_v=2 if k_v=0;
- b_v=1 if k_v=1;
- b_v=0 if k_v>=2.

Find a minimum-cardinality ordinary edge set C satisfying deg_C(v)>=b_v for every vertex. This is polynomial: its complement is a maximum-cardinality b-matching with vertex capacities deg_H(v)-b_v. Negative capacity means infeasibility. General graph capacitated b-matching has a polynomial algorithm. No connectivity constraint is included here.

Set S_0=C union K and L=|C|+|K|. The canonicalization lemma shows L<=OPT on these remaining components. Every vertex has at least two incident resources. Its only defective components, by Section 1, are ordinary simple cycles. Let c be their number. Since each has at least three distinct ordinary edges, c<=|C|/3.

While an ordinary-cycle component remains, add one unused ordinary edge touching it. Such an edge exists by the same ambient-component argument as in Section 2. If it is a chord, it makes the component bicyclic. If it joins another component, it either joins two unicycles, or a unicycle and a robust component, and the joined component is robust. No addition creates a new defective component. Therefore at most c additions suffice.

The returned robust set has size at most

    L+c <= L+|C|/3 <= (4/3)L <= (4/3)OPT.

Combining these solutions with the exact preprocessing solutions preserves 4/3. If an independent fixed charge a>=0 is necessary, the bound a+(4/3)OPT <= (4/3)(a+OPT) holds.

## 4. Strongly connected control realization and path obstruction

Given a connected simple H with at least one edge and no isolated vertices, introduce states p_v for each job, q_e for each ordinary resource, and c_e for each edge. Put arcs c_e -> p_u,p_v,q_e for e=uv and arcs from every p or q state to every c state. This graph is strongly connected. All centers form a free summand of the state transversal matroid, hence loops in its dual. On the p/q states, private q_e columns show that the dual rank of selected Q union J is |J|+nu(Q,V minus J), namely the assignment rank of the ordinary resources Q and one singleton resource per selected p_v. In particular its rank is |V|.

Thus every robust assignment maps directly to a robust actuator placement of identical cardinality, and center actuators are useless for rank. The single source-SCC quota is automatic because every robust spanning set has size at least |V|+1>=3. Consequently the control optimum equals the assignment optimum. The preceding theorem is a 4/3 approximation for this whole control subclass.

When H is a path, retain exactly its two leaf singleton resources. The degree relaxation forces every ordinary path edge. Its output is the path plus its two endpoint singletons, of cardinality |V|+1, which is optimal. Thus it solves exactly the family on which arbitrary basis-plus-exact-repair approaches tend to 3/2.

## 5. Four-source-state hardness family, including source accessibility

For completeness, the control construction is as follows. Given a simple ordinary graph H and distinct jobs s,t, add singleton resources a_s,a_t to its ordinary edge resources. For each resource e introduce states r_e,c_e, and for each job v introduce j_v. Add source states A,B,F_A,F_B. The arcs are r_e -> c_e; c_e -> r_e and c_e -> j_v for each job served by e; the directed cycle A -> F_A -> B -> F_B -> A; the arcs F_B -> r_(a_s), F_A -> r_(a_t); and A -> c_e for every resource. The four source states form the unique source SCC and reach every other state.

The state-transversal columns c_e,F_A,F_B form a free summand, matched respectively from r_e,A,B, so these states are loops of the dual matroid. On the remaining columns, the rows c_e have neighborhoods {r_e} plus their served jobs, and the two filler rows have neighborhoods {A,r_(a_s)} and {B,r_(a_t)}. Without the filler rows and A,B, the dual rank on selected ordinary/singleton resources Q and selected job states J is |J|+nu(Q,V(H) minus J), the assignment rank with a singleton at each selected job. Adding each filler row and its private source column makes that source column a parallel copy of its endpoint resource in the dual. Therefore A and B represent extra singleton copies, and the control matching condition is exactly the assignment condition. This can also be checked directly by expanding the dual rank formula. The independent verifier uses the displayed arcs and direct state-plus-input matchings rather than this identification.

The actual lower-bound-excess construction has a simple ordinary graph H; arbitrary singleton copies at its endpoints and completion singletons at every job; centers/fillers that are loops of the dual matroid; and four source states A,B,F_A,F_B. A and B correspond to additional singleton copies at s and t, whereas F_A,F_B are matroid loops. All other non-source matroid loops can be discarded.

Enumerate all membership patterns Z subset {A,B,F_A,F_B} of size at least two. For that pattern, include exactly Z among these four source states. Interpret selected A/B as mandatory singleton resources and unselected A/B as forbidden singleton resources; their other physical singleton copies remain separate available resources. The selected fillers impose the fixed charge a=|Z intersect {F_A,F_B}|. Run the theorem on the remaining assignment problem, and add the selected filler states. There are only eleven patterns (not sixteen feasible patterns), so this remains polynomial.

For each feasible pattern, the output is at most 4/3 times the optimum restricted to that pattern. An unrestricted optimum has one of these patterns. Taking the smallest returned control placement therefore gives a 4/3 approximation for the full four-source-state control family, without assuming A and B are forced in unrestricted optima. This is an approximation upper bound on the hardness family, not a contradiction to its NP-completeness at d+1 and not a lower bound for the general problem.

## Relationship to the stronger path-cover bound

The companion [path-cover proof](bicircular-pathcover-equivalence.md) establishes an 8/7 approximation for the singleton-completed subclass, including two mandatory singleton resources and the full four-source-state control construction. That bound is stronger than 4/3 on their common scope. The distinctive scope of the present theorem is arbitrary mandatory and forbidden singleton resources; its full polynomial solver is implemented here. The imported Berman-Karpinski path-cover algorithm underlying 8/7 is not implemented in this directory.

## Limitations and next obstruction

Parallel ordinary resources create 2-cycle components and weaken the elementary charging to 3/2. General resource hyperedges serving three or more jobs do not have this degree-cover characterization. Arbitrary source quotas that simultaneously constrain many resource groups are not handled; the constant four-state source gadget is handled by explicit enumeration. Thus extending this approach beyond simple graph-resource assignment requires a new argument, not merely the same counting proof.

## 6. Self-contained polynomial reduction for the degree cover

The degree-cover step does not require trusting an unimplemented b-matching oracle. For capacities c_v=deg_H(v)-b_v, make c_v slot vertices for each original vertex v. For each original edge e=uv make two private port vertices e_u,e_v; join the ports to each other, join e_u to every slot of u, and join e_v to every slot of v. Find a maximum ordinary matching N in this gadget.

An original edge belongs to the excluded capacitated b-matching X exactly when both its ports are matched externally to their endpoint slots. Each vertex participates in at most c_v such excluded edges. Every edge gadget contributes at least one matching edge in a maximum matching: otherwise its internal port edge could be added. A gadget contributes two precisely when both ports are externally matched, and one otherwise. Thus |N|=|E|+|X|. Conversely any capacitated b-matching X gives a gadget matching of size |E|+|X| by assigning its endpoint incidences to distinct slots and internally matching all other edge gadgets. Therefore extracted X is maximum and E minus X is a minimum degree cover.

There are at most 4|E| gadget vertices and O(|E|^2) gadget edges. The existing audited Edmonds blossom routine yields a polynomial implementation, O((|V|+|E|)^3) graph-operation time. `polynomial_degree_cover.py` implements precisely this reduction and copies the existing independently audited blossom function unchanged. The proof above is independent of implementation checks.

## 7. Reproducible audit and nonclaims

Run with standard Python 3.10 or later:

    python verify_four_thirds.py
    python verify_control_four_thirds.py

The first verifier tests all ordinary graphs through four vertices with one optional singleton per vertex and every possible mandatory-loop pattern, plus fixed-seed small cases with zero, one, or two singleton resources per vertex. It independently checks robustness by bipartite job matching rather than the structural characterization. It verifies the exact optimal cost is preserved by the canonical singleton set, compares the actual polynomial blossom-gadget degree cover with exhaustive minimum degree covers and exhaustive maximum capacitated b-matchings, and checks the final 4/3 bound. Pathological cases include isolated vertices, pure cycles, leaves, disconnected graphs, missing singletons, and parallel singleton copies.

The second verifier constructs four-source-state control digraphs, rather than only assignment instances. It evaluates every one of the eleven source membership patterns, checks each algorithm output using full state-plus-input matching and direct accessibility after every actuator failure, and compares its cost with the exact assignment optimum for the same source pattern. It checks arbitrary placements including center states against the independent assignment projection. On the tiny state instances it also directly searches actual control placements up to the exact optimum.

Recorded successful runs:

- 1,896 mandatory-loop fixtures and 94,264 selected-resource masks independently matched and compared with the characterization;
- 18 actual control graphs, 198 source-pattern algorithm outputs, 3,365 direct state/input matching and reachability checks, and 26,203 robust assignment supports in exact-optimum enumerations.

Finite verification supplements the universal proofs; it is not a proof of the ratio or publication novelty. No approximation optimality is established. The ratio for arbitrary one-failure dedicated-actuator structural controllability is still bounded here only by the separate general 3/2 result. This theorem improves the upper bound on the full simple graph-resource control family, including the current exact-hardness reduction and the old repair-tight family, rather than the general problem.

`robust_assignment_solver.py` is the non-enumerative polynomial implementation of the full graph-resource algorithm; `polynomial_degree_cover.py` supplies its blossom gadget. The exhaustive verifier compares this production solver with the independently enumerated optimum on every fixture.
