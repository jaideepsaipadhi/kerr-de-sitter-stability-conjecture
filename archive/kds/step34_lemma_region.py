"""Step 34: the 4-parameter event lemma on the physical region B >= 4t - 1/2 (from c1 >= 4t - 1/2).  Dense scan."""
import numpy as np
from multiprocessing import Pool
from step32_four_param import check

params = [(d, p, b) for d in (0.5, 1.0, 0.25, 2.0) for p in (0.5, 0.2, 0.05, 0.01) for b in (1.0, 0.25, 4.0, 25.0)]

def job(pt):
    t, zf, A, Bfac = pt
    z0 = 1 + zf*(t - 1); B = (4*t - 0.5)*Bfac
    zz = np.linspace(1, z0, 200)
    if np.min(A/zz + B - 3*zz) <= 0:
        return pt, 'Lneg'
    for pr in params:
        ok, info = check(t, z0, A, B, *pr)
        if ok: return pt, 'ok'
    return pt, info

if __name__ == '__main__':
    ts = [1.01, 1.05, 1.2, 1.6, 2.5, 4.0, 7.0, 15.0, 40.0, 100.0, 300.0]
    zfs = [0.001, 0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99, 0.999]
    As = [1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6]
    Bfs = [1.0, 1.3, 3.0, 30.0]
    pts = [(t, zf, A, B) for t in ts for zf in zfs for A in As for B in Bfs]
    with Pool() as pool:
        res = pool.map(job, pts)
    ok = sum(r == 'ok' for _, r in res); neg = sum(r == 'Lneg' for _, r in res)
    fails = [(p, r) for p, r in res if r not in ('ok', 'Lneg')]
    print(f'{ok}/{len(pts) - neg} closed on B >= 4t - 1/2 ({neg} with L <= 0 skipped; none expected)')
    for f in fails[:30]: print('  fail', f)
