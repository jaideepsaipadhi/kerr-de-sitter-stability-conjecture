import numpy as np
import step52d_lift_test as L
res, C4 = L.res, L.C4
ncol = res['ncol']; nkM = res['nk']; compsM = res['comps']; monsKM = res['monsK']; monsM = res['mons']
keepM = [k for k, mo in enumerate(monsM) if sum(mo) <= L.K - 2]; nbM = len(keepM); nb4 = len(L.keep4)
mons4 = L.mons4; keep4 = L.keep4
s0 = 0.8
def eqcol4(ab, mono4): return L.comps4.index(ab)*nb4 + keep4.index(L.idx4[mono4])
def eqcolM(ab, monoM): return compsM.index(ab)*nbM + keepM.index([k for k, mo in enumerate(monsM) if mo == monoM][0])
bad = 0; tot = 0
for pq in compsM:
    for mono in [mo for mo in [monsM[mi] for mi in monsKM]]:
        j = [jj for jj, mi in enumerate(monsKM) if monsM[mi] == mono][0]; row = compsM.index(pq)*nkM + j
        u = [0.0]*(2*ncol); u[ncol + row] = 1.0      # imaginary unit jet
        vre, vim = L.lift(u); a4 = vre @ C4
        for ab in [(0, 0), (1, 1), (3, 3), (0, 1), (0, 3), (1, 3), (0, 2), (1, 2), (2, 2), (2, 3)]:
            # pure r-direction jets: 4-D monomial (0, k, 0, 0) <-> 2-D (k, 0)
            for k in range(0, L.K - 1):
                mv = -float(res['CI'][row, eqcolM(ab, (k, 0))])*(-s0)**sum(1 for i_ in ab if i_ == 2); v4 = a4[eqcol4(ab, (0, k, 0, 0))]
                tot += 1
                if abs(mv - v4) > 1e-8*max(1, abs(mv), abs(v4)):
                    bad += 1
                    if bad <= 8: print(f'MISMATCH h{pq} H-mono{mono} Ric{ab} r^{k}: mode {mv:+.6e} 4D {v4:+.6e}')
print(f'r-direction Ricci jet entries compared: {tot}, mismatches: {bad}')
