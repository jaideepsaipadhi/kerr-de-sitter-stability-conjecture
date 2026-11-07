"""Exact closure counts (K=4, all four Teukolsky conditions) at random rational geometries, points and modes."""
import random, time, numpy as np
from fractions import Fraction as Fr
import step52_exact_modes as S
from step9_freeq import valid_kds
random.seed(7); out = []; t0 = time.time()
while len(out) < 30 and time.time() - t0 < 540:
    Lam = Fr(random.randint(1, 170), 1000); a = Fr(random.randint(5, 140), 100)
    if not valid_kds(float(Lam), float(a)): continue
    rts = np.sort(np.roots([-float(Lam)/3, 0, 1 - float(Lam*a*a)/3, -2, float(a*a)]).real)
    r1, r2 = rts[2], rts[3]
    r0 = Fr(round(r1 + (r2 - r1)*random.uniform(0.1, 0.9), 2)).limit_denominator(100)
    x0 = Fr(random.randint(-9, 9), 10)
    m = random.randint(-4, 4)
    kind = random.choice(['real', 'imag', 'cplx', 'zero'])
    sr = Fr(random.randint(1, 30), 10) if kind in ('real', 'cplx') else Fr(0)
    si = Fr(random.randint(1, 20), 10) if kind in ('imag', 'cplx') else Fr(0)
    if sr == 0 and si == 0 and m == 0: continue          # PD sector checked separately (needs K >= 5)
    res = S.run(Lam, a, Fr(1), r0, x0, (sr, si), m, 4)
    d = S.counts(res, 4)
    out.append((float(Lam), float(a), float(r0), float(x0), (float(sr), float(si)), m, d))
    print(out[-1], flush=True)
print('ALL d_0 = d_1 = 0:', all(o[-1][0] == 0 and o[-1][1] == 0 for o in out), ' samples:', len(out))
