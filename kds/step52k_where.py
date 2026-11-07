import numpy as np
import step52d_lift_test as L
res, C4 = L.res, L.C4; nb4 = len(L.keep4)
best = (0, None)
for bi, b in enumerate(L.basis):
    for part, v in enumerate(L.lift(b)):
        if np.abs(v).max() == 0: continue
        rr_ = v @ C4
        for i in range(10):
            for jj, kk in enumerate(L.keep4):
                val = abs(rr_[i*nb4 + jj])/(np.abs(v).max()*np.abs(C4).max())
                if val > best[0]: best = (val, (bi, part, L.comps4[i], L.mons4[kk], rr_[i*nb4 + jj], np.abs(v).max()))
print('largest Ricci violation:', best)
bi = best[1][0]; b = np.array(L.basis[bi]); ncol = res['ncol']
nz = np.nonzero(b)[0]
print('that basis vector: nonzero entries', len(nz), ' max |entry|', np.abs(b).max())
