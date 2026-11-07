"""Corner limits of the event lemma with the FIXED R3 multiplier (delta, beta, c) = (1, 2, 1/2), on sequences inside R3."""
import numpy as np
from step32_four_param import check
def run(t, zf, A, B):
    x0 = zf*(t - 1); z0 = 1 + x0; L1 = A + B - 3
    ok, info = check(t, z0, A, B, 1.0, 0.5*min(0.9, 2*x0*L1/(t - 1)), 2.0)
    return ok, (round(max(info), 4) if not isinstance(info, str) else info)
fl = lambda t: 5.5*t - 1.25
seqs = {
 'B -> inf':        [(2.0, 0.5, 5.0, fl(2)*10**k) for k in range(0, 8)],
 'A,B -> inf':      [(2.0, 0.5, 10**k, fl(2)*10**k) for k in range(1, 8)],
 'A -> inf':        [(2.0, 0.5, 10**k, fl(2)) for k in range(1, 8)],
 't -> inf':        [(10**k, 0.5, 10**k, fl(10**k)) for k in range(1, 7)],
 'zf -> 0 (in R3)': [(2.0, 10**-k, 10**(2*k + 1), fl(2)) for k in range(1, 6)],
 'zf -> 1':         [(2.0, 1 - 10**-k, 10.0, fl(2)) for k in range(1, 7)],
 't -> 1':          [(1 + 10**-k, 0.5, 10.0, fl(1)) for k in range(1, 7)],
}
for k, seq in seqs.items():
    print(k, [run(*p) for p in seq])
