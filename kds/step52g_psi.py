import numpy as np
import step52d_lift_test as L
res, C4 = L.res, L.C4
ncol = res['ncol']; nkM = res['nk']; compsM = res['comps']; monsKM = res['monsK']; monsM = res['mons']
keepM = [k for k, mo in enumerate(monsM) if sum(mo) <= L.K - 2]; nbM = len(keepM); nb4 = len(L.keep4)
i0m = 10*nbM; i0_4 = 10*nb4          # psi0 block start (mode: complex column block; 4D: re block then im block)
for pq in [(0, 0), (1, 1), (2, 2), (3, 3), (0, 3), (1, 2)]:
    for mono in [(0, 0), (0, 1), (1, 0)]:
        j = [jj for jj, mi in enumerate(monsKM) if monsM[mi] == mono][0]; row = compsM.index(pq)*nkM + j
        u = [0.0]*(2*ncol); u[row] = 1.0
        vre, vim = L.lift(u); a_ = vre @ C4; b_ = vim @ C4
        mv = complex(float(res['CR'][row, i0m]), float(res['CI'][row, i0m]))
        v4c = complex(a_[i0_4], a_[i0_4 + nb4]) + 1j*complex(b_[i0_4], b_[i0_4 + nb4])
        rat = mv/v4c if abs(v4c) > 1e-14 else None
        print(f'h{pq} mono{mono}: mode psi0(p)={mv:.4e}  4D psi0(p)={v4c:.4e}  ratio={rat}')
