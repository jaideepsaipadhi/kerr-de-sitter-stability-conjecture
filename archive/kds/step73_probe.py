import random, time, numpy as np
from fractions import Fraction as Fr
from step72_fast_cert import cert
from step9_freeq import valid_kds
random.seed(11); rows = []
cases = [(Fr(1, 20), Fr(0), Fr(3), Fr(3, 5)), (Fr(1, 20), Fr(1, 1000), Fr(3), Fr(3, 5))]
while len(cases) < 14:
    Lam = Fr(random.randint(1, 175), 1000); a = Fr(random.randint(1, 140), 100)
    if not valid_kds(float(Lam), float(a)): continue
    rts = np.sort(np.roots([-float(Lam)/3, 0, 1 - float(Lam*a*a)/3, -2, float(a*a)]).real)
    r0 = Fr(round(rts[2] + (rts[3] - rts[2])*random.uniform(0.2, 0.8), 2)).limit_denominator(100)
    cases.append((Lam, a, r0, Fr(random.randint(-9, 9), 10)))
for Lam, a, r0, x0 in cases:
    try: d, tb, ta = cert(Lam, a, r0, x0)
    except Exception as e: d = f'ERR {type(e).__name__}'
    print(f'Lam={float(Lam):.3f} a={float(a):.3f} r={float(r0):.2f} x={float(x0):+.1f}: d_4 <= {d}', flush=True)
