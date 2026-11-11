"""Step 42: corner behaviour of the 4-parameter event lemma.  For each corner, follow a sequence into it with a
fixed multiplier family and report the normalised margins (max C/|.|, max D/|.|); a corner is benign when the
margins converge to a strictly negative limit (then the limit problem is the corner lemma)."""
import numpy as np
from step32_four_param import check
MULTS = [(d, p, b) for d in (0.1, 0.25, 0.5, 1.0, 2.0) for p in (0.5, 0.2, 0.05) for b in (0.05, 0.25, 1.0, 4.0, 25.0)]
def best(t, zf, A, B):
    z0 = 1 + zf*(t - 1); out = None
    for mlt in MULTS:
        ok, info = check(t, z0, A, B, *mlt)
        if isinstance(info, str): return 'L<=0'
        sc = max(info)
        if out is None or sc < out[0]: out = (sc, mlt)
    return out
def Bf(t, f): return f*(5.5*t - 1.25)   # B at f times the floor
corners = {
 'A,B -> inf (large m)': [(2.0, 0.5, 10**k, Bf(2, 10**k)) for k in range(0, 6)],
 't -> 1':               [(1 + 10**-k, 0.5, 1.0, Bf(1, 2)) for k in range(1, 7)],
 't -> inf':             [(10**k, 0.5, 1.0*10**k, Bf(10**k, 1.5)) for k in range(1, 6)],
 'zf -> 0':              [(2.0, 10**-k, 1.0, Bf(2, 1.5)) for k in range(1, 7)],
 'zf -> 1':              [(2.0, 1 - 10**-k, 1.0, Bf(2, 1.5)) for k in range(1, 7)],
 'A -> 0':               [(2.0, 0.5, 10**-k, Bf(2, 1.5)) for k in range(1, 7)],
}
for name, seq in corners.items():
    print(name)
    for (t, zf, A, B) in seq:
        print(f'   t={t:<10.6g} zf={zf:<10.6g} A={A:<9.3g} B={B:<9.3g} ->', best(t, zf, A, B), flush=True)
