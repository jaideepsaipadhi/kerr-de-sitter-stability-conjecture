"""zf -> 0 corner with the adaptive p = c * 2 (z0-1) L(1)/(t-1).  Also scan t and A along the corner."""
import numpy as np
from step32_four_param import check
def run(t, zf, A, B):
    z0 = 1 + zf*(t - 1); L1 = A + B - 3; best = None
    for c in (0.05, 0.2, 0.5, 0.9):
        p = c*min(0.9, 2*(z0 - 1)*L1/(t - 1))
        for d in (0.1, 0.5, 2.0):
            for b in (0.05, 1.0, 25.0, 400.0):
                ok, info = check(t, z0, A, B, d, p, b)
                if isinstance(info, str): return info
                if best is None or max(info) < best[0]: best = (max(info), (d, round(p, 8), b))
    return best
for t in (1.1, 2.0, 30.0, 1e3):
    for A in (1e-3, 1.0, 1e3):
        B = 1.5*(5.5*t - 1.25)
        print(f't={t:<6g} A={A:<6g}', [ (zf, round(run(t, zf, A, B)[0], 4)) for zf in (1e-2, 1e-4, 1e-6, 1e-8)], flush=True)
