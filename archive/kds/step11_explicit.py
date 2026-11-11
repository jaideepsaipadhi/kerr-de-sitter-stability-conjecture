"""
Step 11: explicit multipliers, no LP.

Take qt LINEAR on [1,t], through the horizon values qt(1) < 0 < qt(t), with slope q1 > 0.
Then L = (qt^2 - P4)/mu, and the seam is z0 = 1 + (t-1) qt(1)/(qt(1) - qt(t)).

Cosmological side, f+ = 1.  Then (f q')' = 0, so the seam condition holds automatically, and
  C+ = mu' q - 2 mu q1 <= 0 on [z0,t]
  D+ = L' q <= 0 on [z0,t], i.e. L' <= 0 there (q > 0)
  boundary: L(t) > 0 (q(t) > 0)
Event side, f- = (z-1)^p with 0 < p < 1.  Then f(1) = 0, |f'| = o(1/(z-1)), and
(f q')'(z0) = q1 p (z0-1)^{p-1} > 0, so the seam condition holds automatically.
  C- = f [ mu' q - 2 mu q1 - mu q p/(z-1) ] <= 0 on [1,z0]
  D- = L q f' + L' q f + (1/2) q1 (mu f')' <= 0 on [1,z0]
We test p on a grid and report the conditions that fail.

Usage: python3 step11_explicit.py
"""
import numpy as np
from step6_hintz_lp import dual_hintz_form
from step5_window import ctdc_window
from step9_freeq import valid_kds

P = np.polynomial.polynomial


def explicit_check(Lam, a, m, lam, frac, s=-2, M=1.0, ps=(0.25, 0.5, 0.75)):
    lo, hi, _, _ = ctdc_window(Lam, a, M)
    base = dual_hintz_form(Lam, a, M, s, m*(lo + frac*(hi - lo)), m, lam, sign0=-1)
    t, p4, mu = base['t'], base['p4'], base['mu']
    qa, qb = P.polyval(1.0, base['q']), P.polyval(t, base['q'])
    q1 = (qb - qa)/(t - 1)
    q = np.array([qa - q1, q1])
    z0 = 1 + (t - 1)*qa/(qa - qb)
    bump = P.polymul([-1.0, 1.0], [t, -1.0])
    Nq, _ = P.polydiv(P.polysub(P.polymul(q, q), p4), bump)
    L = lambda z: P.polyval(z, Nq)/z
    Lp = lambda z: (P.polyval(z, P.polyder(Nq))*z - P.polyval(z, Nq))/z**2
    mv = lambda z: P.polyval(z, mu); mp_ = lambda z: P.polyval(z, P.polyder(mu)); mpp = lambda z: P.polyval(z, P.polyder(mu, 2))
    qv = lambda z: P.polyval(z, q)
    res = dict(signs_ok=(qa < 0 < qb), t=t, z0=z0)
    # cosmological side
    zc = np.linspace(z0, t, 3001)[1:-1]
    Cp = mp_(zc)*qv(zc) - 2*mv(zc)*q1
    res['cosmo_C'] = float(np.max(Cp/(np.abs(mp_(zc)*qv(zc)) + np.abs(2*mv(zc)*q1))))
    res['cosmo_Lp'] = float(np.max(Lp(zc))/np.max(np.abs(Lp(zc))))
    res['L_t'] = float(L(t - 1e-9))
    # event side
    ze = np.linspace(1, z0, 3001)[1:-1]
    best = None
    for p in ps:
        f = (ze - 1)**p; f1 = p*(ze - 1)**(p - 1); f2 = p*(p - 1)*(ze - 1)**(p - 2)
        Cm = mp_(ze)*f*qv(ze) - 2*mv(ze)*f*q1 - mv(ze)*f1*qv(ze)
        nC = np.abs(mp_(ze)*f*qv(ze)) + np.abs(2*mv(ze)*f*q1) + np.abs(mv(ze)*f1*qv(ze))
        Dm = L(ze)*qv(ze)*f1 + Lp(ze)*qv(ze)*f + 0.5*q1*(mp_(ze)*f1 + mv(ze)*f2)
        nD = np.abs(L(ze)*qv(ze)*f1) + np.abs(Lp(ze)*qv(ze)*f) + 0.5*abs(q1)*(np.abs(mp_(ze)*f1) + np.abs(mv(ze)*f2))
        worst = max(np.max(Cm/nC), np.max(Dm/nD))
        if best is None or worst < best[1]:
            best = (p, worst, float(np.max(Cm/nC)), float(np.max(Dm/nD)))
    res['event_best_p'], res['event_worst'], res['event_C'], res['event_D'] = best
    res['cosmo_ok'] = res['cosmo_C'] <= 0 and res['cosmo_Lp'] <= 0 and res['L_t'] > 0
    res['event_ok'] = best[1] <= 0
    return res


if __name__ == '__main__':
    geoms = [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95), (0.05, 0.9877), (0.08, 0.5), (0.01, 0.99)]
    tot = {'all': 0, 'cosmo': 0, 'event': 0, 'both': 0, 'signs': 0}
    fails = []
    for (Lam, a) in geoms:
        if not valid_kds(Lam, a): print('invalid', Lam, a); continue
        for m in (1, 2, 4, 8):
            for lam in (2.0, 6.0, 20.0, 80.0):
                for fr in (0.02, 0.2, 0.4, 0.6, 0.8, 0.98):
                    r = explicit_check(Lam, a, m, lam, fr)
                    tot['all'] += 1; tot['signs'] += r['signs_ok']
                    tot['cosmo'] += r['cosmo_ok']; tot['event'] += r['event_ok']; tot['both'] += r['cosmo_ok'] and r['event_ok']
                    if not (r['cosmo_ok'] and r['event_ok']):
                        fails.append(((Lam, a, m, lam, fr), {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}))
    print(tot)
    for f in fails[:40]:
        print(f)
