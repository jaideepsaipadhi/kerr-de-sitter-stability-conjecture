"""
Step 18: shape of L2 at ell_min.  Report where L2' has the 'wrong' sign (L2' < 0 on (1,z0) or
L2' > 0 on (z0,t)), and how large the ratio is there compared with the allowed bound H/(mu |z-z0|).
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from step13_largem import leading
from step5_window import ctdc_window
from step17_ellmin import ell_min
from step9_freeq import valid_kds
P = np.polynomial.polynomial

for (Lam, a) in [(0.05, 0.6), (0.02, 0.9), (0.001, 0.95), (0.05, 0.9877), (0.01, 0.99), (0.0001, 0.99), (0.001, 0.999)]:
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    print(f'Lam={Lam} a={a}')
    for wf in (0.02, 0.2, 0.5, 0.8, 0.98):
        lm = ell_min(Lam, a, lo + wf*(hi - lo))
        d = leading(Lam, a, wf, lm)
        t, qa, qb, p4, mu = d['t'], d['qa'], d['qb'], d['p4'], d['mu']
        q1 = (qb - qa)/(t - 1); q = np.array([qa - q1, q1]); z0 = 1 + (t - 1)*qa/(qa - qb)
        bump = P.polymul([-1.0, 1.0], [t, -1.0])
        Nq, _ = P.polydiv(P.polysub(P.polymul(q, q), p4), bump)
        L2 = lambda z: P.polyval(z, Nq)/z
        L2p = lambda z: (P.polyval(z, P.polyder(Nq))*z - P.polyval(z, Nq))/z**2
        muv = lambda z: P.polyval(z, mu); mup = lambda z: P.polyval(z, P.polyder(mu))
        H = lambda z: 2*muv(z) - (z - z0)*mup(z)
        ze = np.linspace(1, z0, 4001)[1:-1]; zc = np.linspace(z0, t, 4001)[1:-1]
        ratio_e = (-L2p(ze)/L2(ze))/(H(ze)/(muv(ze)*(z0 - ze)))   # need <= 1
        ratio_c = (L2p(zc)/L2(zc))/(H(zc)/(muv(zc)*(zc - z0)))     # need <= 1
        zmax = np.linspace(1, t, 20001)[np.argmax(L2(np.linspace(1, t, 20001)))]
        print(f'   wfrac={wf}: z0={z0:.4f} argmax L2={zmax:.4f} | event: frac wrong-sign {np.mean(L2p(ze) < 0):.2f}, '
              f'max ratio {ratio_e.max():+.3f} | cosmo: frac wrong-sign {np.mean(L2p(zc) > 0):.2f}, max ratio {ratio_c.max():+.3f}')
