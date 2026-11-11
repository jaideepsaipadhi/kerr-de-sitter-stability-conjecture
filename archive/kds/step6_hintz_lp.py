"""
Step 6: Hintz's two-current argument on the real dual, posed as LP feasibility.

The dual in Hintz form is
    (-d mu d + Lt(z)) R -/+ ... = 0  with  mu W = qt^2 - Lt mu,
where qt is a real quadratic and Lt is linear.  The boundary conditions fix the signs of
qt(1) and qt(t).  Inside the CTdC window these signs are opposite, so qt has a zero z0
in (1,t).  The leftover freedom is the sign of qt(0).

Identity (rederived here for a general quadratic q and linear L):
  J' = C^f |u'|^2 + D^f |u|^2,
  C^f = mu' f q - 2 mu f q' - mu f' q,
  D^f = L q f' + L' q f + (1/2)(mu (f q')')',
  J   = -mu f q |u'|^2 - mu f q' Re(u' conj u) + ((1/2) mu (f q')' + L f q)|u|^2.
At a horizon J -> L f q |B|^2.  At the seam, J(z0) = -mu f q' X + (1/2) mu (f q')' |R|^2,
with X = Re(R' conj R)(z0).

Conditions that give R = 0 (both orientations of q'(z0) are covered):
  cosmological side [z0,t]: f > 0, C <= 0, D <= 0, (f q')'(z0) <= 0, and L(t) q(t) > 0
  event side [1,z0]: f(1) = 0, f >= 0, |f'| = o(1/(z-1)), C <= 0, D <= 0, (f q')'(z0) >= 0
For a fixed sign of q'(z0), the cosmological side gives X of one sign (strictly, unless
B_c = 0) and the event side gives X of the other sign.  So X = 0, hence B_c = 0, hence R = 0.

Each side becomes an LP over polynomial f = sum c_k phi_k on a sample grid.  The objective
maximises a slack tau (normalised).  tau >= 0 means the side is feasible on the grid.
"""
import itertools
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.optimize import linprog
from step3_dual import teukolsky_heun, heun_params_from_operator, euler_dual, normal_form_Q, real_branch
from step5_window import ctdc_window

P = np.polynomial.polynomial


def dual_hintz_form(Lam, a, M, s, omega, m, lam, sign0=+1):
    data = teukolsky_heun(Lam, a, M, s, complex(omega), m, lam)
    t = float(complex(data['t']).real)
    e = data['exps']
    h = heun_params_from_operator(data['Pz'], data['Qz'], data['z'], data['t'],
                                  e['0'][0], e['1'][0], e['t'][0],
                                  {'e': 2*s + 1, 'zinf': sp.Float(str(data['zinf']), 35)}, s)
    _, ht = real_branch(h, complex(t))
    A = ht['A']
    mu = np.array([0.0, -t, 1 + t, -1.0])  # z(z-1)(t-z) ascending coefficients
    mup = P.polyder(mu); mupp = P.polyder(mup)
    def P4(z):
        Q = normal_form_Q(ht, complex(t), complex(z)).real
        return P.polyval(z, mu)**2*Q - (P.polyval(z, mup)**2 - 2*P.polyval(z, mu)*P.polyval(z, mupp))/4
    zs = np.linspace(1.02, t - 0.02, 11)
    p4 = P.polyfit(zs, [P4(z) for z in zs], 4)
    q1 = (1j*(t - 1)*(A - h['d'])/2).real
    qtt = (1j*abs(t - t**2)*(h['e'] - A)/2).real
    P40 = P.polyval(0.0, p4)
    assert P40 >= -1e-9, 'P4(0) < 0: no real quadratic qt'
    # check consistency: qt(1)^2 = P4(1), qt(t)^2 = P4(t)
    c1 = abs(q1**2 - P.polyval(1.0, p4)); ct = abs(qtt**2 - P.polyval(t, p4))
    q = P.polyfit([0.0, 1.0, t], [sign0*np.sqrt(max(P40, 0)), q1, qtt], 2)
    rem = P.polysub(P.polymul(q, q), p4)  # = Lt mu
    Lt, r_ = P.polydiv(rem, mu)
    return dict(t=t, mu=mu, q=q, L=Lt, p4=p4, consistency=(c1, ct, float(np.max(np.abs(r_)))))


def side_lp(form, z0, side, K=10, N=300):
    t, mu, q, L = form['t'], form['mu'], form['q'], form['L']
    a_, b_ = (z0, t) if side == 'cosmo' else (1.0, z0)
    # basis: shifted Legendre on [a_, b_]; for the event side multiply by (z-1) so f(1) = 0
    x = lambda z: 2*(z - a_)/(b_ - a_) - 1
    Leg = np.polynomial.legendre.Legendre
    zm1 = Leg.fit([a_, b_], [a_ - 1.0, b_ - 1.0], 1, domain=[a_, b_])  # the polynomial z - 1, native basis
    def basis_poly(k):
        c = np.zeros(k + 1); c[k] = 1
        Lk = Leg(c, domain=[a_, b_])
        return Lk*zm1 if side == 'event' else Lk
    def basis(k, z, der):
        p_ = basis_poly(k)
        return p_.deriv(der)(z) if der else p_(z)
    zz = np.linspace(a_, b_, N)[1:-1]
    qv, qp, qpp = P.polyval(zz, q), P.polyval(zz, P.polyder(q)), P.polyval(zz, P.polyder(q, 2))
    muv, mup = P.polyval(zz, mu), P.polyval(zz, P.polyder(mu))
    Lv, Lp = P.polyval(zz, L), P.polyval(zz, P.polyder(L))
    rows_C, rows_D, rows_f = [], [], []
    for k in range(K):
        f0, f1, f2 = basis(k, zz, 0), basis(k, zz, 1), basis(k, zz, 2)
        C = mup*f0*qv - 2*muv*f0*qp - muv*f1*qv
        fq1 = f1*qp + f0*qpp
        fq2 = f2*qp + 2*f1*qpp
        D = Lv*qv*f1 + Lp*qv*f0 + 0.5*(mup*fq1 + muv*fq2)
        rows_C.append(C); rows_D.append(D); rows_f.append(f0)
    RC, RD, RF = np.array(rows_C).T, np.array(rows_D).T, np.array(rows_f).T
    sC = np.max(np.abs(RC)) or 1; sD = np.max(np.abs(RD)) or 1
    # seam derivative (f q')'(z0)
    seam = np.array([basis(k, z0, 1)*P.polyval(z0, P.polyder(q)) + basis(k, z0, 0)*P.polyval(z0, P.polyder(q, 2)) for k in range(K)])
    # variables: c_0..c_{K-1}, tau ; maximise tau
    A_ub, b_ub = [], []
    for i in range(len(zz)):
        A_ub.append(list(RC[i]/sC) + [1.0]); b_ub.append(0.0)  # C/sC + tau <= 0
        A_ub.append(list(RD[i]/sD) + [1.0]); b_ub.append(0.0)  # D/sD + tau <= 0
        A_ub.append(list(-RF[i]) + [1.0]); b_ub.append(0.0)    # f >= tau
    if side == 'cosmo':
        A_ub.append(list(seam) + [0.0]); b_ub.append(0.0)      # (f q')'(z0) <= 0
        A_eq = [list(np.array([basis(k, t, 0) for k in range(K)])) + [0.0]]  # normalise f(t) = 1
    else:
        A_ub.append(list(-seam) + [0.0]); b_ub.append(0.0)     # (f q')'(z0) >= 0
        A_eq = [list(np.array([basis(k, z0, 0) for k in range(K)])) + [0.0]]  # f(z0) = 1
    b_eq = [1.0]
    bounds = [(-1e3, 1e3)]*K + [(None, 1.0)]
    res = linprog(c=[0]*K + [-1.0], A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
    return (res.x[-1] if res.success else None), res


def test_frequency(Lam, a, M, s, omega, m, lam, K=10):
    out = []
    for sign0 in (+1, -1):
        form = dual_hintz_form(Lam, a, M, s, omega, m, lam, sign0)
        t = form['t']
        zr = [r.real for r in P.polyroots(form['q']) if abs(r.imag) < 1e-12 and 1 < r.real < t]
        if len(zr) != 1:
            out.append(dict(sign0=sign0, note=f'{len(zr)} zeros of qt in (1,t)'))
            continue
        z0 = zr[0]
        Lq_t = P.polyval(t, form['L'])*P.polyval(t, form['q'])
        tc, _ = side_lp(form, z0, 'cosmo', K)
        te, _ = side_lp(form, z0, 'event', K)
        out.append(dict(sign0=sign0, z0=z0, Lq_at_t=Lq_t, tau_cosmo=tc, tau_event=te,
                        Lcoef=form['L'], qcoef=form['q'], consistency=form['consistency']))
    return out


if __name__ == '__main__':
    s = -2
    for (Lam, a, M) in [(0.05, 0.6, 1.0), (0.02, 0.9, 1.0)]:
        lo, hi, Om2, Om1 = ctdc_window(Lam, a, M)
        print(f'\n##### Lam={Lam} a={a} M={M}   CTdC window omega/m in ({lo:.5f},{hi:.5f})')
        for m in (1, 2, 3):
            for lam in (2.0, 6.0, 20.0):
                for frac in (0.2, 0.5, 0.8):
                    om = m*(lo + frac*(hi - lo))
                    for r in test_frequency(Lam, a, M, s, om, m, lam):
                        if 'note' in r:
                            print(f'm={m} lam={lam:5.1f} frac={frac}: sign0={r["sign0"]:+d}  {r["note"]}')
                        else:
                            tc = 'infeas' if r['tau_cosmo'] is None else f"{r['tau_cosmo']:+.2e}"
                            te = 'infeas' if r['tau_event'] is None else f"{r['tau_event']:+.2e}"
                            print(f"m={m} lam={lam:5.1f} frac={frac}: sign0={r['sign0']:+d} z0={r['z0']:.4f} "
                                  f"L(t)q(t)={r['Lq_at_t']:+.3e} tau_cosmo={tc} tau_event={te}  "
                                  f"L={np.round(r['Lcoef'],3)} cons={max(r['consistency']):.1e}")
