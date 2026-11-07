"""
Step 9: let qt vary.

The two-current identity only needs
  * qt smooth and real on [1,t], with the fixed horizon values qt(1), qt(t) (these set the
    exponents, i.e. the boundary conditions), and
  * L := (qt^2 - P4)/mu, which is then smooth on [1,t] because qt(z_h)^2 = P4(z_h) at z_h = 1, t.
qt does NOT have to be the quadratic interpolant through z = 0, and L does not have to be linear.

Family:  qt = q_lin + (z-1)(t-z) g(z), where q_lin is the line through the two horizon values
and g is a polynomial.  For each g we run the piecewise LP of step 8 on both sides, and
maximise min(tau_cosmo, tau_event) over g with Nelder-Mead.
A frequency counts as closed only if the best (qt, f+-) pair passes the fine-grid check
and L(t) qt(t) > 0.

Usage:  python3 step9_freeq.py quick    # 7 hardest cases
        python3 step9_freeq.py medium   # 168 cases, parallel; ~20-40 min on an 8-core laptop
        python3 step9_freeq.py full     # 672 cases; several hours, only if medium is clean
"""
import sys
import numpy as np
from scipy.optimize import minimize
from step6_hintz_lp import dual_hintz_form
from step5_window import ctdc_window
import step8_piecewise as S8

P = np.polynomial.polynomial


def make_form(base, gcoef):
    t, p4 = base['t'], base['p4']
    q1 = P.polyval(1.0, base['q']); qt_ = P.polyval(t, base['q'])
    qlin = np.array([q1 - (qt_ - q1)/(t - 1), (qt_ - q1)/(t - 1)])
    bump = P.polymul([-1.0, 1.0], [t, -1.0])              # (z-1)(t-z)
    q = P.polyadd(qlin, P.polymul(bump, gcoef))
    num = P.polysub(P.polymul(q, q), p4)
    Nq, rem = P.polydiv(num, bump)                        # exact division up to rounding
    def Lf(z): return P.polyval(z, Nq)/z
    def Lpf(z): return (P.polyval(z, P.polyder(Nq))*z - P.polyval(z, Nq))/z**2
    return dict(t=t, mu=base['mu'], q=q, Lf=Lf, Lpf=Lpf, p4=p4, divrem=float(np.max(np.abs(rem))))


def coeff_arrays(form, zz):
    q = form['q']
    return dict(q=P.polyval(zz, q), qp=P.polyval(zz, P.polyder(q)), qpp=P.polyval(zz, P.polyder(q, 2)),
                mu=P.polyval(zz, form['mu']), mup=P.polyval(zz, P.polyder(form['mu'])),
                L=form['Lf'](zz), Lp=form['Lpf'](zz))


S8.coeff_arrays = coeff_arrays  # reuse step 8 machinery with a general L


def zero_in(form):
    t = form['t']
    zz = np.linspace(1, t, 4001)
    v = P.polyval(zz, form['q'])
    idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    if len(idx) != 1:
        return None
    from scipy.optimize import brentq
    return brentq(lambda z: P.polyval(z, form['q']), zz[idx[0]], zz[idx[0] + 1])


def score(base, g, cfg=(8, 4, 1.0)):
    form = make_form(base, np.asarray(g))
    z0 = zero_in(form)
    if z0 is None:
        return -1.0, None
    t = form['t']
    Lq = form['Lf'](t - 1e-9)*P.polyval(t, form['q'])
    if Lq <= 0:
        return -1.0, None
    tc, solc = S8.piecewise_lp(form, z0, 'cosmo', *cfg[:2], grading=cfg[2])
    te, sole = S8.piecewise_lp(form, z0, 'event', *cfg[:2], grading=cfg[2])
    tc = -1.0 if tc is None else tc
    te = -1.0 if te is None else te
    return min(tc, te), dict(form=form, z0=z0, tc=tc, te=te, solc=solc, sole=sole, Lq=Lq)


def valid_kds(Lam, a, M=1.0):
    import numpy as _np
    rts = _np.roots([-Lam/3, 0, 1 - Lam*a**2/3, -2*M, a**2])
    if _np.max(_np.abs(rts.imag)) > 1e-10*max(1, _np.max(_np.abs(rts))):
        return False
    r3, r0, r1, r2 = _np.sort(rts.real)
    # subextremal KdS = four distinct real roots r3 < 0 < r0 < r1 < r2 (NO Lam M^2 < 1/9 condition: that is the a=0 bound only)
    return r3 < 0 < r0 < r1 - 1e-6 and r1 < r2 - 1e-3


def run(Lam, a, m, lam, frac, s=-2, M=1.0, gdeg=2):
    if not valid_kds(Lam, a, M):
        return 'INVALID', dict(best=float('nan'))
    lo, hi, _, _ = ctdc_window(Lam, a, M)
    base = dual_hintz_form(Lam, a, M, s, m*(lo + frac*(hi - lo)), m, lam, sign0=-1)
    # start from g that reproduces the quadratic interpolant (deg 0 correction), then optimise
    t = base['t']
    q1 = P.polyval(1.0, base['q']); qt_ = P.polyval(t, base['q'])
    qlin = np.array([q1 - (qt_ - q1)/(t - 1), (qt_ - q1)/(t - 1)])
    bump = P.polymul([-1.0, 1.0], [t, -1.0])
    g0, _ = P.polydiv(P.polysub(base['q'], qlin), bump)
    g0 = np.concatenate([np.atleast_1d(g0), np.zeros(gdeg + 1)])[:gdeg + 1]
    best = (score(base, g0), g0)
    starts = [g0] + [g0 + np.random.default_rng(k).normal(0, 0.5*max(1, np.abs(g0).max()), gdeg + 1) for k in range(4)]
    for st in starts:
        r = minimize(lambda g: -score(base, g)[0], st, method='Nelder-Mead',
                     options=dict(maxfev=120, xatol=1e-4, fatol=1e-6))
        sc = score(base, r.x)
        if sc[0] > best[0][0]:
            best = (sc, r.x)
        if best[0][0] > 1e-3:
            break
    (val, info), g = best
    if info is None or val < 0:
        return 'OPEN', dict(best=val)
    fcs = {}
    for side, key in (('cosmo', 'solc'), ('event', 'sole')):
        sol = info[key]; extra = {}
        for it in range(8):  # cutting planes
            fc, bad = S8.fine_check(info['form'], info['z0'], side, sol, return_points=True)
            if not bad:
                break
            for i_, pts in bad.items():
                extra[i_] = np.concatenate([extra.get(i_, np.array([])), pts])
            tau, sol2 = S8.piecewise_lp(info['form'], info['z0'], side, 8, 4, grading=1.0, extra=extra)
            if tau is None or tau < 0:
                break
            sol = sol2
        fcs[side] = fc
    fcC, fcE = fcs['cosmo'], fcs['event']
    ok = all(fc['worstC'] <= 1e-10 and fc['worstD'] <= 1e-10 and fc['minf_rel'] >= -1e-10 for fc in (fcC, fcE))
    return ('CLOSED' if ok else 'OPEN'), dict(best=val, g=g, z0=info['z0'], fineC=fcC, fineE=fcE, divrem=info['form']['divrem'])


def _job(case):
    try:
        st, info = run(*case)
    except Exception as ex:
        st, info = 'ERROR', dict(best=float('nan'), err=str(ex))
    return case, st, info.get('best', float('nan'))


if __name__ == '__main__':
    import multiprocessing as mp_
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    geoms = [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95), (0.05, 0.9877), (0.08, 0.5), (0.01, 0.99)]
    if mode == 'quick':
        cases = [(0.05, 0.6, 1, 2.0, 0.05), (0.05, 0.6, 1, 2.0, 0.95), (0.02, 0.9, 1, 50.0, 0.95),
                 (0.02, 0.9, 1, 2.0, 0.95), (0.001, 0.95, 2, 10.0, 0.35), (0.001, 0.95, 2, 10.0, 0.65),
                 (0.001, 0.95, 1, 2.0, 0.95)]
    elif mode == 'medium':
        cases = [(L_, a_, m, lam, fr) for (L_, a_) in geoms for m in (1, 4) for lam in (2.0, 50.0) for fr in (0.03, 0.5, 0.97)]
    else:
        cases = [(L_, a_, m, lam, fr) for (L_, a_) in geoms for m in (1, 2, 4, 8)
                 for lam in (2.0, 6.0, 20.0, 80.0) for fr in (0.02, 0.2, 0.4, 0.6, 0.8, 0.98)]
    with mp_.Pool() as pool:
        n = c = 0
        for case, st, best in pool.imap_unordered(_job, cases):
            n += st != 'INVALID'; c += st == 'CLOSED'
            print(f"{st:7} {case}  best min-tau={best:+.2e}", flush=True)
    print(f'\nSUMMARY: {c}/{n} closed (valid geometries only)')
