# A 5/3 obstruction for maximum-coverage joint repair

## Maximum-coverage greedy can approach 5/3

For d=3k, introduce targets p_0,...,p_{d-1}, q_0,...,q_{d-2} and centers c_0,...,c_{d-2}. Edges are c_j→p_j,p_{j+1},q_j, and every target→every center. The graph has n=3d−2 vertices and is strongly connected. At most d−1 targets and d−1 centers can be matched, and these bounds are attained simultaneously, so its matching deficiency is d.

The set P={p_i} is a size-d spanning set and satisfies the source quota. Its repair supports are exactly A(q_j)={p_j,p_{j+1}}; centers cover nothing. Indeed matching all q except q_j and adding p_i is possible precisely for i=j or j+1, since each q has its unique center predecessor.

W={p_0,p_{d-1}}∪{q_0,...,q_{d-2}} is robust and has size d+1. If an endpoint fails, the d−1 unactuated path vertices match one-to-one along the path. If q_j fails, match q_j from c_j; match internal p_1,...,p_j from c_0,...,c_{j-1}, and p_{j+1},...,p_{d-2} from c_{j+1},...,c_{d-2}. Match all center targets independently from distinct target left copies. These are state matchings after every actuator failure. Strong connectivity supplies accessibility. The rank lower bound proves OPT=d+1.

Give q_1,q_4,...,q_{3k-2} the lowest labels among q. Maximum-coverage greedy chooses these k pair covers first. Remaining uncovered points are p_0,p_3,...,p_{3k-3}; every remaining cover hits at most one. Another k repairs are necessary. Its output has d+2k=5d/3 points, and ratio 5d/[3(d+1)] tends to 5/3. Thus this greedy variant does not have a 3/2 guarantee, despite all small random tests.

The executable numbering is q first (with the indicated tie order), then centers, then p. The original PR's deterministic max-flow solver returns exactly P on all tested family members, so the phenomenon does not require changing its relaxation tie-breaking. `python audit_path_family.py` independently constructs the graph, checks strong connectivity, verifies explicit matching and robust-optimum certificates, and invokes `joint_solver.py` for d=3,6,...,39. It checks the actual flow placement, target repair supports, output size and residual size without importing a matching oracle. In particular d=12 gives 20 actuators versus OPT=13, so the implemented deterministic greedy is not a 3/2 approximation. The asymptotic statement allows the indicated optimum relaxation and tie choices; deterministic flow behavior has been checked on the listed finite range.

Optimal fixed-P repair on this family is a minimum edge cover of a d-vertex path, of size ceil(d/2), so the exact joint-repair ratio tends to 3/2. The 5/3 obstruction is specific to greedy choice, not to joint optimization, and establishes neither a 7/4 lower bound nor hardness for the original problem.
