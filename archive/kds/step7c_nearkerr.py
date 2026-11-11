"""Near-Kerr probe: is the cosmological-side failure an LP artefact?  Look at Hintz's
analogue f+ = 1/|q'| directly, then rerun the LP with a larger basis and wider bounds."""
import numpy as np
from step6_hintz_lp import dual_hintz_form, side_lp
from step5_window import ctdc_window
P = np.polynomial.polynomial
for Lam, a in [(0.001, 0.95), (0.05, 0.6)]:
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    for frac in (0.35, 0.65):
        m, lam = 2, 10.0
        form = dual_hintz_form(Lam, a, 1.0, -2, m*(lo + frac*(hi - lo)), m, lam, sign0=-1)
        t, mu, q, L = form['t'], form['mu'], form['q'], form['L']
        z0 = [r.real for r in P.polyroots(q) if abs(r.imag) < 1e-12 and 1 < r.real < t][0]
        zz = np.linspace(z0, t, 4001)[1:-1]
        qv, qp = P.polyval(zz, q), P.polyval(zz, P.polyder(q)); qpp = P.polyval(zz, P.polyder(q, 2))
        muv, mup = P.polyval(zz, mu), P.polyval(zz, P.polyder(mu))
        Lv, Lp = P.polyval(zz, L), P.polyval(zz, P.polyder(L))
        f = 1/np.abs(qp); f1 = -np.sign(qp)*qpp/qp**2
        C = mup*f*qv - 2*muv*f*qp - muv*f1*qv
        D = Lv*qv*f1 + Lp*qv*f   # (f q')' = 0 for this f, so the mu-term drops out
        print(f'Lam={Lam} frac={frac}: t={t:.4f} z0={z0:.4f}  q coeffs={np.round(q,4)}  L coeffs={np.round(L,3)}')
        print(f'    f+=1/|q\'|: max C={C.max():+.3e} (min {C.min():+.3e}), max D={D.max():+.3e}; sign changes of q\' on [z0,t]: {int(np.sum(np.diff(np.sign(qp))!=0))}')
        print(f'    L(z) on [z0,t]: min {Lv.min():+.3f} max {Lv.max():+.3f};  q(t)={P.polyval(t,q):+.4f}')
        for K in (12, 20):
            tau, res = side_lp(form, z0, 'cosmo', K=K, N=800)
            print(f'    LP cosmo K={K}: tau={tau}')
