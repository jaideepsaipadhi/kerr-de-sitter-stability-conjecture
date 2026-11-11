"""
Step 7: stress test of the LP result.

For each test frequency (sign0 = -1):
  1. Solve the LP on a grid of N = 300 points, as in step 6.
  2. Re-evaluate C^f, D^f and f for the returned multiplier on a 6000-point grid that
     reaches within 1e-6 of each endpoint, and report the worst normalised violation.
  3. Check the endpoint data: L(t) q(t) > 0, the sign of (f q')'(z0), f(1) = 0 and the
     behaviour of f' near z = 1.

Geometries: generic, near-extremal rotation (a close to its maximum for the given Lambda),
large Lambda, and small Lambda (close to Kerr).

Usage:  python3 step7_verify.py quick   (~5 min)
        python3 step7_verify.py full    (9 geometries x 5 m x 5 lambda x 13 freqs; roughly 1-2 h)
"""
import sys
import numpy as np
from step6_hintz_lp import dual_hintz_form, side_lp
from step5_window import ctdc_window
from step3_dual import kds_roots

P = np.polynomial.polynomial


def eval_side(form, z0, side, res, K, N=6000):
    t, mu, q, L = form['t'], form['mu'], form['q'], form['L']
    a_, b_ = (z0, t) if side == 'cosmo' else (1.0, z0)
    c = res.x[:K]
    Leg = np.polynomial.legendre.Legendre
    zm1 = Leg.fit([a_, b_], [a_ - 1.0, b_ - 1.0], 1, domain=[a_, b_])
    f = Leg([0.0], domain=[a_, b_])
    for k in range(K):
        e = np.zeros(k + 1); e[k] = 1
        Lk = Leg(e, domain=[a_, b_])
        if side == 'event':
            Lk = Lk*zm1
        f = f + c[k]*Lk
    eps = 1e-6*(b_ - a_)
    zz = np.concatenate([np.linspace(a_ + eps, a_ + 0.01*(b_ - a_), N//4),
                         np.linspace(a_ + 0.01*(b_ - a_), b_ - 0.01*(b_ - a_), N//2),
                         np.linspace(b_ - 0.01*(b_ - a_), b_ - eps, N//4)])
    f0, f1, f2 = f(zz), f.deriv(1)(zz), f.deriv(2)(zz)
    qv, qp, qpp = P.polyval(zz, q), P.polyval(zz, P.polyder(q)), P.polyval(zz, P.polyder(q, 2))
    muv, mup = P.polyval(zz, mu), P.polyval(zz, P.polyder(mu))
    Lv, Lp = P.polyval(zz, L), P.polyval(zz, P.polyder(L))
    C = mup*f0*qv - 2*muv*f0*qp - muv*f1*qv
    D = Lv*qv*f1 + Lp*qv*f0 + 0.5*(mup*(f1*qp + f0*qpp) + muv*(f2*qp + 2*f1*qpp))
    sC, sD = np.max(np.abs(C)), np.max(np.abs(D))
    seam = f.deriv(1)(z0)*P.polyval(z0, P.polyder(q)) + f(z0)*P.polyval(z0, P.polyder(q, 2))
    seam_scale = abs(f.deriv(1)(z0)*P.polyval(z0, P.polyder(q))) + abs(f(z0)*P.polyval(z0, P.polyder(q, 2)))
    return dict(seam_scale=seam_scale, maxC=np.max(C)/sC, maxD=np.max(D)/sD, minf=np.min(f0)/np.max(np.abs(f0)),
                seam=seam, f_at_1=f(1.0), fprime_at_1=f.deriv(1)(1.0))


def run_case(Lam, a, M, m, lam, frac, s=-2, K=12):
    lo, hi, _, _ = ctdc_window(Lam, a, M)
    om = m*(lo + frac*(hi - lo))
    form = dual_hintz_form(Lam, a, M, s, om, m, lam, sign0=-1)
    t = form['t']
    zr = [r.real for r in P.polyroots(form['q']) if abs(r.imag) < 1e-12 and 1 < r.real < t]
    if len(zr) != 1:
        return dict(status=f'{len(zr)} zeros of qt in (1,t)')
    z0 = zr[0]
    out = dict(z0=z0, Lq_t=P.polyval(t, form['L'])*P.polyval(t, form['q']), cons=max(form['consistency']))
    for side in ('cosmo', 'event'):
        tau, res = side_lp(form, z0, side, K)
        if tau is None or tau < 0:
            out[side] = dict(tau=tau)
            continue
        v = eval_side(form, z0, side, res, K)
        v['tau'] = tau
        out[side] = v
    ok = (out['Lq_t'] > 0)
    for side, sgn in (('cosmo', -1), ('event', +1)):
        v = out[side]
        if v.get('tau') is None or v['tau'] < 0 or 'maxC' not in v:
            ok = False; continue
        ok &= v['maxC'] <= 1e-9 and v['maxD'] <= 1e-9 and v['minf'] >= -1e-9 and sgn*v['seam'] >= -1e-7*(1 + abs(v['seam_scale']))
    out['status'] = 'CLOSED' if ok else 'OPEN'
    return out


def max_a(Lam, M=1.0):
    # largest a with four real roots (subextremal), by bisection
    lo, hi = 0.0, 3.0
    for _ in range(60):
        mid = (lo + hi)/2
        try:
            rts = kds_roots(Lam, mid, M)
            good = len(rts) == 4 and rts[1] > 0 and rts[1] < rts[2] - 1e-9
        except Exception:
            good = False
        lo, hi = (mid, hi) if good else (lo, mid)
    return lo


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    geoms = [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95)]
    amax = max_a(0.05)
    geoms.append((0.05, 0.97*amax))
    if mode == 'full':
        geoms += [(0.08, 0.5), (0.1, 0.2), (0.01, 0.99*max_a(0.01)), (0.1, 0.95*max_a(0.1))]
        ms, lams, fracs = (1, 2, 3, 5, 8), (2.0, 4.0, 10.0, 40.0, 120.0), np.linspace(0.02, 0.98, 13)
    else:
        ms, lams, fracs = (1, 2, 4), (2.0, 10.0, 50.0), (0.05, 0.35, 0.65, 0.95)
    tot = closed = 0
    for (Lam, a) in geoms:
        try:
            rts = kds_roots(Lam, a, 1.0)
        except Exception as ex:
            print('skip', Lam, a, ex); continue
        print(f'\n### Lambda={Lam}  a={a:.5f}  roots={[round(float(x), 4) for x in rts]}')
        for m in ms:
            for lam in lams:
                for frac in fracs:
                    try:
                        r = run_case(Lam, a, 1.0, m, lam, frac)
                    except Exception as ex:
                        r = dict(status=f'ERROR {ex}')
                    tot += 1; closed += r['status'] == 'CLOSED'
                    if r['status'] != 'CLOSED':
                        print(f'  m={m} lam={lam} frac={frac:.2f}: {r["status"]}  ', {k: v for k, v in r.items() if k != 'status'})
        print(f'  running total: {closed}/{tot} closed')
    print(f'\nSUMMARY: {closed}/{tot} test frequencies closed by the two-current argument on the dual')
