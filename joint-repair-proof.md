# Joint repair gives a 7/4 approximation

## Scope
Unit costs, at most one dedicated input per state, arbitrary directed state graph, and robustness against one input failure. Retain the minimum quota-spanning relaxation and notation from [the factor-two note](two-approximation.md): M=T*, rank d, minimum relaxation P of size L, and D the coloops of M|P. Source SCCs have at least two vertices; otherwise the instance is infeasible.

## Repair sets
For b in D define H_b=cl_M(P minus {b}) and C_b=E minus H_b. Since P spans M and b is a coloop of M|P, H_b is a rank-(d-1) hyperplane and C_b a cocircuit. Moreover, C_b intersects P in precisely {b}. For e outside P define A_e={b in D: e in C_b}. Equivalently, b belongs to A_e exactly when P minus {b} plus {e} spans M.

Each b belongs to at least one A_e because M has no coloops. For any R outside P, P union R is one-failure spanning if and only if the A_e for e in R cover D. Necessity: if C_b is missed by R, deleting b leaves the whole placement inside H_b. Sufficiency: deleting a repair retains P; deleting b in D retains an individual replacement spanning set; deleting a noncritical original element retains P minus that element. Source quotas are preserved by addition.

## Algorithm
Set U=D and R empty. While some e outside P satisfies |A_e intersection U|>=2, add e to R and remove A_e from U. Any such choice is valid; selecting maximum coverage is a deterministic natural choice with vertex-index tie breaking. Once no such e exists, write t=|U| and select an individual repair for each remaining b in U.

If q repairs were chosen in the first phase, they covered at least 2q distinct elements of D; thus 2q<=|D|-t. The last phase uses t distinct repairs, since no external element covers two residual elements. Therefore |R|=q+t<=(|D|+t)/2.

## Packing lower bound
For distinct b,c in residual U, C_b and C_c are disjoint. Their intersections with P are the distinct singletons {b},{c}. An element outside P in their intersection would have |A_e intersection U|>=2, contradicting the stopping condition.

Every one-failure spanning solution O intersects every cocircuit C in at least two elements. If O misses C then O does not span; if it meets C in a single element, deletion of that element leaves a set contained in the complementary hyperplane. Consequently the t pairwise disjoint residual cocircuits give OPT>=2t. This bound applies to arbitrary robust solutions, not only supersets of P.

## Guarantee
The returned solution S satisfies

|S|=L+|R| <= L+(|D|+t)/2 <= L+d/2+t/2 <= (7/4)OPT.

The final inequality uses L<=OPT, d<=OPT, and 2t<=OPT. If d>0, the stronger rank bound d+1<=OPT yields |S|<=(7/4)OPT-1/2. If d=0 then D is empty and the quota relaxation is already exact.

A computable lower bound is max(L,d+1,2t) for d>0, and L for d=0. Do not confuse this with an exact optimum or with a proof of 3/2. The general additive certificate is |S|<=L+(|D|+t)/2.

## Cubic implementation
Use the matching F covering I=E minus P provided by the relaxation flow. For each b in D, perform alternating reachability starting at the unmatched right b, going right-to-left along a state edge and left-to-right along the matching edge. No unmatched left vertex is reachable because b is critical. An external right vertex e belongs to the reached set if and only if swapping the exposed endpoint from b to e yields a matching covering I minus {e} plus {b}, equivalently b in A_e. Thus all A_e can be computed by at most n searches of O(n^2), without repeated rank or matching calls. The equivalence follows either from the alternating path flip, or from the alternating component in the symmetric difference between F and a matching covering I-e+b.

The greedy stage makes at most n selections and scans at most n sets of at most n elements per selection, hence O(n^3). Together with the baseline flow, this is O(n^3) time and O(n^2) space in the graph RAM model. The reference implementation uses Python sets; the bound counts standard graph/set operations and does not claim bit-cost independence for arbitrary-precision Python integer masks.

## Verification
Use Python 3.10 or newer. `joint_solver.py` reuses the existing `two_approx_solver.py` and `bounded_core_solver.py`; `check_joint.py` uses `verify_core_theorem.py`. Run `python check_joint.py`. It keeps the relaxation matching certificate, original critical elements, coverage incidence, residual packing, and chosen repair states. `check_joint.py` exhaustively tested all 66,066 directed graphs with self-loops on 1 through 4 vertices plus 600 fixed-seed random graphs on 5 through 8 vertices. All 66,666 passed. There were 38,951 feasible instances. It independently checks optimum, robustness, packing lower bound, repair-count bound, and 7/4. Worst observed ratio 1.5 is experimental only.

## Literature and limitations
Poças, Pequito and Ramos, Automatica 190 (2026), 113074, https://doi.org/10.1016/j.automatica.2026.113074, treats the same no-duplicate dedicated model as its primary Problem 2 and states factors 5/2 for general digraphs and 2 for strongly connected digraphs. Bentert, Fomin, Golovach and Morelle, Fault-Tolerant Matroid Bases, ESA 2025, https://arxiv.org/html/2506.22010v1, treats the abstract robust-spanning problem and its parameterized complexity. Neither source establishes the present quota-first joint-repair bound in the passages inspected. This limited search does not establish publication novelty.

The result does not establish a 3/2 factor, weighted or multiple-failure guarantees, or tightness of 7/4. No re-selection or exchange of original P is needed. Maximum coverage is not required for the proof; arbitrary coverage of at least two residual critical elements suffices.

### Important prior 7/4 analysis
Adjiashvili, Bindewald and Michaels, *Robust Assignments with Vulnerable Nodes*, https://arxiv.org/pdf/1703.06074, Theorem 7 (pp. 16–18 in the PDF), already establishes a 1.75 approximation for uniform cardinality vertex-robust assignment. Its proof distinguishes greedy steps covering at least two failures from single-cover steps and derives a disjoint-neighborhood lower bound for the latter. The numerical constant and productive/nonproductive greedy analysis must therefore not be presented as new. The present argument adapts that framework to the dual-transversal matroid with source-SCC quotas, using global cocircuits in place of job neighborhoods and a quota-first spanning relaxation in place of a matching basis. The cited theorem is not directly applicable simply by identifying a dual-transversal matroid with a transversal one. Whether this extension or the same controllability guarantee appears elsewhere remains unestablished.

## Independent audit commands
Run `python audit_seven_fourths.py` to check all optimum relaxation sets, residual cocircuit disjointness, packing lower bounds and the guarantee, using subset-DP rank calculations independent of the solver matching. Run `python audit_joint_implementation.py` to check the actual solver, full repair coverage, exact relaxation, residual disjointness, robust output and iterator inputs. Each audit covers all 66,066 looped digraphs through four vertices plus 800 fixed-seed random graphs on five through eight vertices (seed 20261104). These finite checks supplement the proof.
