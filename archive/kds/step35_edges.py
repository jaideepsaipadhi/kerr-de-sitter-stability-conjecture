"""Step 35: rerun the step-34 failures with p scaled to the near-horizon constraint p = c (z0-1) L(1)/mu'(1) and a wider (delta, beta) set."""
import numpy as np
from multiprocessing import Pool
from step32_four_param import check
from step34_lemma_region import job as job34

def job(pt):
    t, zf, A, Bfac = pt
    z0 = 1 + zf*(t - 1); B = (4*t - 0.5)*Bfac
    L1 = A + B - 3; mup1 = t - 1
    pst = 2*(z0 - 1)*L1/mup1
    for c in (0.5, 0.2, 0.05):
        p = min(0.9, c*pst)
        for d in (0.1, 0.25, 0.5, 1.0, 2.0, 4.0):
            for b in (0.05, 0.25, 1.0, 4.0, 25.0, 100.0):
                ok, info = check(t, z0, A, B, d, p, b)
                if ok: return pt, ('ok', c, d, b)
    return pt, info

if __name__ == '__main__':
    ts = [1.01, 1.05, 1.2, 1.6, 2.5, 4.0, 7.0, 15.0, 40.0, 100.0, 300.0]
    zfs = [0.001, 0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99, 0.999]
    As = [1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6]
    Bfs = [1.0, 1.3, 3.0, 30.0]
    pts = [(t, zf, A, B) for t in ts for zf in zfs for A in As for B in Bfs]
    with Pool() as pool:
        first = pool.map(job34, pts)
        fails = [p for p, r in first if r != 'ok']
        second = pool.map(job, fails)
    still = [(p, r) for p, r in second if r[0] != 'ok']
    print(f'after the adaptive retry: {len(pts) - len(still)}/{len(pts)} closed; {len(still)} remaining')
    for s in still[:20]: print('  fail', s)
