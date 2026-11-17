"""Explicit strongly connected family tight for independent repair choices."""
from two_approx_solver import solve
from bounded_core_solver import matching_number, source_components

def instance(d):
    assert d >= 2
    n = 3*d+2; rows = [0]*n
    h,g = d,3*d+1
    for i in range(d):
        e,f = d+1+i,2*d+1+i
        rows[i] = (1<<e)|(1<<f)
        rows[e] = rows[f] = (1<<i)|(1<<h)
    rows[h] = sum(1<<j for j in range(d+1,n))
    rows[g] = sum(1<<j for j in range(d+1))
    return rows

for d in range(2,26):
    rows=instance(d); n=len(rows); full=(1<<n)-1
    assert source_components(rows)==[full]
    assert matching_number(rows,full)==n-d
    witness=tuple(range(d+1,2*d+1))+(3*d+1,)
    mask=sum(1<<i for i in witness)
    for i in witness:
        rem=mask^(1<<i)
        assert matching_number(rows,full^rem)==n-rem.bit_count()
    result=solve(rows)
    assert result['relaxation_value']==d
    assert len(result['actuators'])==2*d
    assert len(result['critical_states'])==d
    assert len(result['repair_states'])==d
print('PASS: d=2,...,25; OPT=d+1 certified by witness and rank lower bound; algorithm returns 2d')
