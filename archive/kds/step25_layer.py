"""
Step 25: explicit finite-m multipliers with the seam boundary layer.  No LP.

Exact finite-m data from the dual (step 6): qt linear through the horizon values, L = (qt^2 - P4)/mu.
Cosmological side: f+ = 1.  Conditions: C+ = -q1 H <= 0 (automatic), D+ = q L' <= 0 (needs L' <= 0), L(t) q(t) > 0.
Event side: f'/f = h(z) = (1+delta) g(z) (1 - exp(-kappa m^4 (z0-z)^2)) + p (1-u)/(z-1),
  with g = -L'/L, u = (z-1)/(z0-1), kappa = (L(z0)/m^2)/mu(z0).
  f is obtained by integrating h; f' = h f and f'' = (h' + h^2) f.
  Check C <= 0 and D <= 0 on a grid refined geometrically toward z0 (it resolves the 1/m^2 layer)
  and toward z = 1.  At the seam h(z0) = 0, so f'(z0) = 0 >= 0 (equality, allowed).

Usage: python3 step25_layer.py
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from step6_hintz_lp import dual_hintz_form
from step5_window import ctdc_window
P = np.polynomial.polynomial


def build(Lam, a, m, ell, wf, s=-2):
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    b = dual_hintz_form(Lam, a, 1.0, s, m*(lo + wf*(hi - lo)), m, ell*m*m, sign0=-1)
    t, p4, mu = b['t'], b['p4'], b['mu']
    qa, qb = P.polyval(1.0, b['q']), P.polyval(t, b['q'])
    q1 = (qb - qa)/(t - 1); q = np.array([qa - q1, q1]); z0 = 1 + (t - 1)*qa/(qa - qb)
    bump = P.polymul([-1.0, 1.0], [t, -1.0])
    N, _ = P.polydiv(P.polysub(P.polymul(q, q), p4), bump)
    L = lambda z: P.polyval(z, N)/z
    Lp = lambda z: (P.polyval(z, P.polyder(N))*z - P.polyval(z, N))/z**2
    Lpp = lambda z: (P.polyval(z, P.polyder(N, 2))*z**2 - 2*P.polyval(z, P.polyder(N))*z + 2*P.polyval(z, N))/z**3
    return dict(t=t, q=q, q1=q1, z0=z0, mu=mu, L=L, Lp=Lp, Lpp=Lpp, m=m)


def check(F, delta=0.5, p=0.5):
    t, z0, q1, m = F['t'], F['z0'], F['q1'], F['m']
    mu = lambda z: P.polyval(z, F['mu']); mup = lambda z: P.polyval(z, P.polyder(F['mu']))
    qv = lambda z: P.polyval(z, F['q'])
    L, Lp, Lpp = F['L'], F['Lp'], F['Lpp']
    out = {}
    # cosmological side, f = 1
    zc = np.linspace(z0, t, 20001)[1:-1]
    out['cosmo_Lp_max'] = float(np.max(Lp(zc)))
    out['Lq_t'] = float(L(t - 1e-12)*qv(t))
    # event side
    kappa = (L(z0)/m**2)/mu(z0)
    x0 = z0 - 1
    sgrid = np.concatenate([np.geomspace(1e-14*x0, 1e-2*x0, 6000), np.linspace(1e-2*x0, x0*(1 - 1e-9), 6000)])
    z = z0 - sgrid
    z = np.sort(np.concatenate([z, 1 + np.geomspace(1e-12*x0, 1e-2*x0, 3000)]))
    z = z[(z > 1) & (z < z0)]
    g = -Lp(z)/L(z); gp = -(Lpp(z)*L(z) - Lp(z)**2)/L(z)**2
    E = np.exp(-kappa*m**4*(z0 - z)**2); Ep = E*(2*kappa*m**4*(z0 - z))
    psi = (1 + delta)*(1 - E); psip = -(1 + delta)*Ep
    u = (z - 1)/x0
    hB = p*(1 - u)/(z - 1); hBp = -p/(z - 1)**2
    h = g*psi + hB; hp = gp*psi + g*psip + hBp
    # C/f and D/f (f > 0 cancels); f'/f = h, f''/f = h' + h^2
    qv_, muv, mupv, Lv, Lpv = qv(z), mu(z), mup(z), L(z), Lp(z)
    C = mupv*qv_ - 2*muv*q1 - muv*h*qv_
    nC = np.abs(mupv*qv_) + np.abs(2*muv*q1) + np.abs(muv*h*qv_)
    D = Lv*qv_*h + Lpv*qv_ + 0.5*q1*(mupv*h + muv*(hp + h**2))
    nD = np.abs(Lv*qv_*h) + np.abs(Lpv*qv_) + 0.5*abs(q1)*(np.abs(mupv*h) + np.abs(muv*(hp + h**2)))
    out['event_C_max'] = float(np.max(C/nC)); out['event_D_max'] = float(np.max(D/nD))
    out['zD'] = float(z[np.argmax(D/nD)]); out['z0'] = z0
    out['ok'] = out['cosmo_Lp_max'] <= 0 and out['Lq_t'] > 0 and out['event_C_max'] <= 0 and out['event_D_max'] <= 0
    return out


if __name__ == '__main__':
    cases = [(0.05, 0.6, 0.5, 1.2), (0.001, 0.95, 0.2, 1.0), (0.01, 0.99, 0.5, 0.9), (0.05, 0.9877, 0.9, 0.8),
             (0.05, 0.6, 0.02, 1.0), (0.05, 0.6, 0.98, 1.0), (0.001, 0.999, 0.2, 0.92)]
    for (Lam, a, wf, ell) in cases:
        line = []
        for m in (8, 32, 128, 512):
            r = check(build(Lam, a, m, ell, wf))
            line.append(f"m={m}: {'OK ' if r['ok'] else 'NO '} C={r['event_C_max']:+.1e} D={r['event_D_max']:+.1e}"
                        f" L'c={r['cosmo_Lp_max']:+.1e}" + ('' if r['ok'] else f" @z={r['zD']:.6f}/z0={r['z0']:.6f}"))
        print(f'Lam={Lam} a={a} wf={wf} ell={ell}:\n    ' + '\n    '.join(line), flush=True)
