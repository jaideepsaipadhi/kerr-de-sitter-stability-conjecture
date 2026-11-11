"""
Step 13: the large-m (semiclassical) limit of the two-current conditions.

Scaling: omega = m w, lambda = m^2 ell, m -> infinity.  Then qt = m q^(1) + O(1),
P4 = m^2 P2 + O(m), L = m^2 L2 + O(m).

Limiting conditions, with qt linear:
  cosmological side (f+ = 1):   mu' q - 2 mu q1 <= 0,   L2' <= 0  on [z0,t],   L2(t) > 0
  event side:  a valid f exists iff
      Phi(z) := (2 mu q1 - mu' q)/(mu |q|) + L2'/L2 >= 0  on (1, z0),  and  L2 > 0 there.
These come from C <= 0 (an upper bound on f'/f) and the O(m^3) part of D <= 0 (a lower
bound on f'/f).

The leading-order objects are computed by evaluating at two large m values and
extrapolating (checking that the extrapolation has converged).

Usage: python3 step13_largem.py
"""
import numpy as np
import mpmath as mp
from step6_hintz_lp import dual_hintz_form
from step5_window import ctdc_window
from step9_freeq import valid_kds

P = np.polynomial.polynomial
mp.mp.dps = 40


def leading(Lam, a, wfrac, ell, s=-2, M=1.0, ms=(400, 800)):
    lo, hi, _, _ = ctdc_window(Lam, a, M)
    w = lo + wfrac*(hi - lo)
    data = []
    for m in ms:
        b = dual_hintz_form(Lam, a, M, s, m*w, m, ell*m*m, sign0=-1)
        t = b['t']
        qa, qb = P.polyval(1.0, b['q'])/m, P.polyval(t, b['q'])/m
        data.append((m, t, qa, qb, np.array(b['p4'])/m**2))
    (m1, t, qa1, qb1, p1), (m2, _, qa2, qb2, p2) = data
    # Richardson-type extrapolation assuming an O(1/m) correction
    qa = 2*qa2 - qa1; qb = 2*qb2 - qb1; p4 = 2*p2 - p1
    conv = max(abs(qa2 - qa1), abs(qb2 - qb1))/max(abs(qa2), abs(qb2))
    return dict(t=t, qa=qa, qb=qb, p4=p4, mu=b['mu'], conv=conv)


def check(Lam, a, wfrac, ell):
    d = leading(Lam, a, wfrac, ell)
    t, qa, qb, p4, mu = d['t'], d['qa'], d['qb'], d['p4'], d['mu']
    out = dict(conv=d['conv'], signs=(qa < 0 < qb))
    if not out['signs']:
        return out
    q1 = (qb - qa)/(t - 1); q = np.array([qa - q1, q1]); z0 = 1 + (t - 1)*qa/(qa - qb)
    bump = P.polymul([-1.0, 1.0], [t, -1.0])
    Nq, _ = P.polydiv(P.polysub(P.polymul(q, q), p4), bump)
    L2 = lambda z: P.polyval(z, Nq)/z
    L2p = lambda z: (P.polyval(z, P.polyder(Nq))*z - P.polyval(z, Nq))/z**2
    muv = lambda z: P.polyval(z, mu); mup = lambda z: P.polyval(z, P.polyder(mu))
    qv = lambda z: P.polyval(z, q)
    zc = np.linspace(z0, t, 4001)[1:-1]
    C = mup(zc)*qv(zc) - 2*muv(zc)*q1
    out['cosmo_C'] = float(np.max(C/(np.abs(mup(zc)*qv(zc)) + 2*np.abs(muv(zc)*q1))))
    out['cosmo_L2p'] = float(np.max(L2p(zc))/np.max(np.abs(L2p(zc))))
    out['L2_t'] = float(L2(t - 1e-9))
    ze = np.linspace(1, z0, 4001)[1:-1]
    up = (2*muv(ze)*q1 - mup(ze)*qv(ze))/(muv(ze)*np.abs(qv(ze)))
    lw = -L2p(ze)/L2(ze)
    out['event_L2min'] = float(np.min(L2(ze)))
    out['event_Phi_rel'] = float(np.min((up - lw)/(np.abs(up) + np.abs(lw))))
    out['cosmo_ok'] = out['cosmo_C'] <= 0 and out['cosmo_L2p'] <= 0 and out['L2_t'] > 0
    out['event_ok'] = out['event_L2min'] > 0 and out['event_Phi_rel'] >= 0
    return out


if __name__ == '__main__':
    geoms = [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95), (0.05, 0.9877), (0.08, 0.5), (0.01, 0.99)]
    for (Lam, a) in geoms:
        if not valid_kds(Lam, a):
            print('invalid', Lam, a); continue
        rows = []; nc = ne = n = 0; worst_conv = 0
        for wf in (0.02, 0.2, 0.4, 0.6, 0.8, 0.98):
            for ell in (1.0, 1.5, 3.0, 10.0, 100.0):
                r = check(Lam, a, wf, ell); n += 1
                worst_conv = max(worst_conv, r['conv'])
                if not r['signs']:
                    rows.append((wf, ell, 'SIGNS')); continue
                nc += r['cosmo_ok']; ne += r['event_ok']
                if not (r['cosmo_ok'] and r['event_ok']):
                    rows.append((wf, ell, {k: round(v, 3) for k, v in r.items() if isinstance(v, float)}))
        print(f'Lam={Lam} a={a}: cosmo {nc}/{n}, event {ne}/{n}, worst m-convergence {worst_conv:.1e}')
        for row in rows[:12]:
            print('    fail', row)
