"""
Step 12: event side, Hintz Case II style.

Take qt linear and L = (qt^2 - P4)/mu as in step 11.  On [1,z0] we need f >= 0, f(1) = 0,
C <= 0, D <= 0 and (f q')'(z0) = q1 f'(z0) >= 0.  Near z = 1 the equation D f = 0 has the
regular exponent p* = 2 L(1)|q(1)| / (q1 mu'(1)).  We try:
  (i)  f = (z-1)^{c p*} for c in (0.3, 0.6, 0.9)
  (ii) the solution of D f = 0 with f ~ (z-1)^{p*} at z = 1 (D = 0 exactly), for which we
       check C <= 0, f > 0 and f'(z0) >= 0
"""
import numpy as np
from scipy.integrate import solve_ivp
from step6_hintz_lp import dual_hintz_form
from step5_window import ctdc_window
from step9_freeq import valid_kds
P = np.polynomial.polynomial

def setup(Lam, a, m, lam, frac, s=-2, M=1.0):
    lo, hi, _, _ = ctdc_window(Lam, a, M)
    base = dual_hintz_form(Lam, a, M, s, m*(lo + frac*(hi - lo)), m, lam, sign0=-1)
    t, p4, mu = base['t'], base['p4'], base['mu']
    qa, qb = P.polyval(1.0, base['q']), P.polyval(t, base['q'])
    q1 = (qb - qa)/(t - 1); q = np.array([qa - q1, q1]); z0 = 1 + (t - 1)*qa/(qa - qb)
    bump = P.polymul([-1.0, 1.0], [t, -1.0])
    Nq, _ = P.polydiv(P.polysub(P.polymul(q, q), p4), bump)
    F = dict(t=t, z0=z0, q1=q1, qa=qa,
             q=lambda z: P.polyval(z, q), mu=lambda z: P.polyval(z, mu),
             mup=lambda z: P.polyval(z, P.polyder(mu)), mupp=lambda z: P.polyval(z, P.polyder(mu, 2)),
             L=lambda z: P.polyval(z, Nq)/z, Lp=lambda z: (P.polyval(z, P.polyder(Nq))*z - P.polyval(z, Nq))/z**2)
    F['pstar'] = 2*F['L'](1.0)*abs(qa)/(q1*F['mup'](1.0))
    return F

def CD(F, z, f, f1, f2):
    q, mu, mup, L, Lp, q1 = F['q'](z), F['mu'](z), F['mup'](z), F['L'](z), F['Lp'](z), F['q1']
    C = mup*f*q - 2*mu*f*q1 - mu*f1*q
    nC = abs(mup*f*q) + abs(2*mu*f*q1) + abs(mu*f1*q) + 1e-300
    D = L*q*f1 + Lp*q*f + 0.5*q1*(mup*f1 + mu*f2)
    nD = abs(L*q*f1) + abs(Lp*q*f) + 0.5*abs(q1)*(abs(mup*f1) + abs(mu*f2)) + 1e-300
    return C/nC, D/nD

def try_power(F, p):
    """f = g(u), u = (z-1)/(z0-1), g(u) = u^p (1 + p - p u): g(0) = 0, g' >= 0, g'(1) = 0."""
    x0 = F['z0'] - 1
    z = np.linspace(1, F['z0'], 6001)[1:-1]; u = (z - 1)/x0
    g = u**p*(1 + p - p*u)
    g1 = p*(1 + p)*u**(p - 1)*(1 - u)
    g2 = p*(1 + p)*((p - 1)*u**(p - 2)*(1 - u) - u**(p - 1))
    C, D = CD(F, z, g, g1/x0, g2/x0**2)
    return max(C.max(), D.max())

def try_ode(F):
    # D f = 0:  (1/2) q1 mu f'' + ((1/2) q1 mu' + L q) f' + L' q f = 0.  Substitute f = x^p* h.
    p = F['pstar']; z0 = F['z0']
    eps = 1e-6*(z0 - 1)
    # Start slightly inside with the Frobenius leading term (adequate for a check).
    def rhs(z, y):
        f, f1 = y
        q, mu, mup, L, Lp, q1 = F['q'](z), F['mu'](z), F['mup'](z), F['L'](z), F['Lp'](z), F['q1']
        f2 = -((0.5*q1*mup + L*q)*f1 + Lp*q*f)/(0.5*q1*mu)
        return [f1, f2]
    x0 = eps
    sol = solve_ivp(rhs, [1 + x0, z0], [x0**p, p*x0**(p - 1)], rtol=1e-10, atol=1e-14, dense_output=True, method='LSODA')
    if not sol.success: return None
    z = np.linspace(1 + x0, z0, 4001)
    f, f1 = sol.sol(z)
    f2 = np.array([rhs(zz, [ff, gg])[1] for zz, ff, gg in zip(z, f, f1)])
    C, _ = CD(F, z, f, f1, f2)
    return dict(maxC=float(C.max()), minf=float(f.min()), fprime_z0=float(f1[-1]))

if __name__ == '__main__':
    geoms = [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95), (0.05, 0.9877), (0.08, 0.5), (0.01, 0.99)]
    tot = ok_pow = ok_ode = ok_any = 0; fails = []; pst = []
    for (Lam, a) in geoms:
        if not valid_kds(Lam, a): continue
        for m in (1, 2, 4, 8):
            for lam in (2.0, 6.0, 20.0, 80.0):
                for fr in (0.02, 0.2, 0.4, 0.6, 0.8, 0.98):
                    F = setup(Lam, a, m, lam, fr); tot += 1; pst.append(F['pstar'])
                    wp = min(try_power(F, p) for p in (0.1, 0.25, 0.5, 0.75, 0.9))
                    o = try_ode(F)
                    po = wp <= 0; oo = o is not None and o['maxC'] <= 1e-10 and o['minf'] > 0 and o['fprime_z0'] >= 0
                    ok_pow += po; ok_ode += oo; ok_any += po or oo
                    if not (po or oo): fails.append(((Lam, a, m, lam, fr), round(F['pstar'], 4), round(wp, 4), o))
    print(f'total {tot}: power OK {ok_pow}, ODE (Case II) OK {ok_ode}, either {ok_any}')
    print(f'p* range: [{min(pst):.4f}, {max(pst):.4f}],  p*<1 in {sum(p < 1 for p in pst)}/{tot}')
    for f in fails[:25]: print(f)
