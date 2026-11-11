"""
Step 14: large-m cosmological side with a general multiplier f+.
A valid f+ exists iff  Psi(z) = -L2'/L2 - (mu' q - 2 mu q1)/(mu q) >= 0 on (z0,t), with L2 > 0 there.
Compare the relative margin with the f+ = 1 margin (max of normalised L2') in the thin cases.
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from step13_largem import leading
from step9_freeq import valid_kds
P = np.polynomial.polynomial

def psi_margin(Lam, a, wf, ell):
    d = leading(Lam, a, wf, ell)
    t, qa, qb, p4, mu = d['t'], d['qa'], d['qb'], d['p4'], d['mu']
    q1 = (qb - qa)/(t - 1); q = np.array([qa - q1, q1]); z0 = 1 + (t - 1)*qa/(qa - qb)
    bump = P.polymul([-1.0, 1.0], [t, -1.0])
    Nq, _ = P.polydiv(P.polysub(P.polymul(q, q), p4), bump)
    z = np.linspace(z0, t, 6001)[1:-1]
    L2 = P.polyval(z, Nq)/z; L2p = (P.polyval(z, P.polyder(Nq))*z - P.polyval(z, Nq))/z**2
    muv, mup, qv = P.polyval(z, mu), P.polyval(z, P.polyder(mu)), P.polyval(z, q)
    up = -L2p/L2; lw = (mup*qv - 2*muv*q1)/(muv*qv)
    rel = (up - lw)/(np.abs(up) + np.abs(lw))
    f1margin = np.max(L2p)/np.max(np.abs(L2p))
    return float(rel.min()), float(np.min(L2)), float(f1margin), float(z[np.argmin(rel)]), z0, t

for (Lam, a) in [(0.001, 0.95), (0.01, 0.99), (0.05, 0.9877), (0.0001, 0.99), (0.001, 0.999)]:
    if not valid_kds(Lam, a): print('invalid', Lam, a); continue
    for wf in (0.5, 0.9, 0.99, 0.999):
        for ell in (1.0, 100.0):
            r, L2min, f1, zw, z0, t = psi_margin(Lam, a, wf, ell)
            print(f'Lam={Lam} a={a} wfrac={wf} ell={ell}: Psi_rel_min={r:+.4f} (at z={zw:.3f}; z0={z0:.3f}, t={t:.2f})  L2min={L2min:+.3e}   [f+=1 margin {f1:+.4f}]')
