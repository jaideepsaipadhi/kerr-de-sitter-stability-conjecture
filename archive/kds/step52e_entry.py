import numpy as np
from fractions import Fraction as Fr
import step52_exact_modes as S
import step52d_lift_test as L          # reuses res (sigma from argv default 0) and 4-D objects
res, C4 = L.res, L.C4
ncol = res['ncol']; nkM = res['nk']; compsM = res['comps']; monsKM = res['monsK']; monsM = res['mons']
keepM = [k for k, mo in enumerate(monsM) if sum(mo) <= L.K - 2]; nbM = len(keepM)
nb4 = len(L.keep4)
for (p, q) in [(1, 1), (0, 0), (0, 3), (3, 3), (0, 1)]:
    u = [0.0]*(2*ncol); u[compsM.index((p, q))*nkM + 0] = 1.0     # H_pq = 1 (degree 0), real
    v4 = L.lift(u)[0]                                             # real part of e^{i phi} H
    for (a_, b_) in [(1, 1), (0, 0), (3, 3)]:
        jm = compsM.index((a_, b_))*nbM                          # degree-0 entry of mode Ric_ab
        mode_val = complex(float(res['CR'][compsM.index((p, q))*nkM, jm]), float(res['CI'][compsM.index((p, q))*nkM, jm]))
        j4 = L.comps4.index((a_, b_))*nb4
        val4 = (v4 @ C4)[j4]
        print(f'h_{p}{q}=1 -> Ric_{a_}{b_}(p):  mode {mode_val.real:+.6e}  4D {val4:+.6e}')
