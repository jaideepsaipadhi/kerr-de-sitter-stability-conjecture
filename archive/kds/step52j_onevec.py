import numpy as np
import step52d_lift_test as L
res, C4 = L.res, L.C4; ncol = res['ncol']
nb4 = len(L.keep4); i_r1 = L.keep4.index(L.idx4[(0, 1, 0, 0)])
for bi in range(3):
    b = np.array(L.basis[bi])
    vre, vim = L.lift(list(b))
    direct = (vre @ C4)[L.comps4.index((1, 1))*nb4 + i_r1]
    # prediction from unit jets: sum_k b_k * lift(e_k)
    pred = 0.0
    for k in np.nonzero(b)[0]:
        e = [0.0]*(2*ncol); e[k] = 1.0
        pred += b[k]*(L.lift(e)[0] @ C4)[L.comps4.index((1, 1))*nb4 + i_r1]
    print(f'basis[{bi}]: max|b| = {np.abs(b).max():.2e},  4D Ric_rr at r^1: direct {direct:+.6e}, sum of unit lifts {pred:+.6e}')
