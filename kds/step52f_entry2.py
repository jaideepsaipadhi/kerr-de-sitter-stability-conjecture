import numpy as np
import step52d_lift_test as L
res, C4 = L.res, L.C4
ncol = res['ncol']; nkM = res['nk']; compsM = res['comps']; monsKM = res['monsK']; monsM = res['mons']
keepM = [k for k, mo in enumerate(monsM) if sum(mo) <= L.K - 2]; nbM = len(keepM); nb4 = len(L.keep4)
def unit(pq, mono):
    u = [0.0]*(2*ncol); j = [jj for jj, mi in enumerate(monsKM) if monsM[mi] == mono][0]
    u[compsM.index(pq)*nkM + j] = 1.0; return u, compsM.index(pq)*nkM + j
for pq in [(2, 2), (1, 2), (0, 2), (2, 3), (1, 1), (0, 0)]:
    for mono in [(0, 0), (0, 1), (1, 0), (0, 2)]:
        u, row = unit(pq, mono)
        v4 = L.lift(u)[0]; r4 = v4 @ C4
        for ab in [(1, 1), (0, 0), (0, 3)]:
            jm = compsM.index(ab)*nbM; j4 = L.comps4.index(ab)*nb4
            mv = float(res['CR'][row, jm])
            if abs(mv - r4[j4]) > 1e-9*max(1, abs(mv)): print(f'MISMATCH h{pq} mono{mono} -> Ric{ab}(p): mode {mv:+.6e}  4D {r4[j4]:+.6e}')
print('done')
