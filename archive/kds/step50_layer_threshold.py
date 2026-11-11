"""Step 50: does the layer multiplier cover all A >= A_R2 (so that R2 + layer = everything with p > 1)?
For each (t, zf) test A = A_R2 * {1, 3, 10, 1e2, 1e4} and B = floor * {1, 10, 1e3}, plus the R1 boundary.
Multiplier: step32 check(t, z0, A, B, delta, p, beta), adaptive p = c*min(0.9, 2 x0 L(1)/(t-1))."""
import numpy as np
from step32_four_param import check
MU = [(d, b, c) for d in (0.05, 0.2, 0.5, 1.0) for b in (0.25, 1.0, 4.0, 25.0) for c in (0.2, 0.5, 0.9)]
res = []; fails = []
for t in np.geomspace(1.01, 3e3, 14):
    for zf in (0.05, 0.2, 0.4, 0.6, 0.8, 0.95):
        x0 = zf*(t - 1); z0 = 1 + x0
        AR2 = max(z0*((t - 1)/(2*x0**2) - 3), 1e-3)
        Bf = 5.5*t - 1.25
        for Af in (1, 3, 10, 1e2, 1e4):
            for Bm in (1, 10, 1e3):
                A, B = AR2*Af, Bf*Bm
                L1 = A + B - 3
                if 2*x0*L1/(t - 1) <= 1: continue          # R1 covers it
                if 2*x0**2*(A/z0 + 3) < t - 1: continue    # R2 covers it
                okk = False
                for (d, b, c) in MU:
                    p = c*min(0.9, 2*x0*L1/(t - 1))
                    ok, info = check(t, z0, A, B, d, p, b)
                    if ok: okk = True; break
                res.append(okk)
                if not okk: fails.append((round(t, 3), zf, f'{A:.3g}', f'{B:.3g}'))
print(f'R3 samples: {len(res)}, layer multiplier closes {sum(res)}')
for f in fails[:20]: print('  fail', f)
