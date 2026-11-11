import numpy as np
from step32_four_param import check
pts = []
for t in np.geomspace(1.01, 3e3, 14):
    for zf in (0.05, 0.2, 0.4, 0.6, 0.8, 0.95):
        x0 = zf*(t - 1); z0 = 1 + x0; AR2 = max(z0*((t - 1)/(2*x0**2) - 3), 1e-3); Bf = 5.5*t - 1.25
        for Af in (1, 3, 10, 1e2, 1e4):
            for Bm in (1, 10, 1e3):
                A, B = AR2*Af, Bf*Bm; L1 = A + B - 3
                if 2*x0*L1/(t - 1) <= 1 or 2*x0**2*(A/z0 + 3) < t - 1: continue
                pts.append((t, z0, A, B, x0, L1))
best = []
for d in (0.05, 0.1, 0.2, 0.5, 1.0):
    for b in (0.25, 0.5, 1.0, 2.0, 4.0):
        for c in (0.2, 0.5, 0.9):
            n = 0; worst = -9
            for (t, z0, A, B, x0, L1) in pts:
                ok, info = check(t, z0, A, B, d, c*min(0.9, 2*x0*L1/(t - 1)), b)
                n += ok
                if not isinstance(info, str): worst = max(worst, max(info))
            best.append((n, -worst, d, b, c))
best.sort(reverse=True)
print(len(pts), 'points; top universal choices (closed, -worst margin, delta, beta, c):')
for r in best[:6]: print('  ', r)
