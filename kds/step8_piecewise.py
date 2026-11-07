"""
Step 8: piecewise-polynomial multipliers for the two-current argument on the real dual.

Setup as in step 6.  Each side's interval is split into n_pieces pieces, and f is a
Legendre polynomial on each piece.
  * f is continuous at the breakpoints (a jump in f puts an indefinite cross term into J).
  * At a breakpoint b, J jumps by (1/2) mu(b) q'(b) [f'](b) |u|^2.  We need this <= 0, so
    J stays non-increasing.  The condition is q'(b) (f'_R - f'_L) <= 0.
Interior conditions are unchanged: f >= tau, C <= -tau, D <= -tau (normalised).  The seam
and endpoint conditions are as in step 6.

After solving, the multiplier is re-checked on a fine grid (2000 points per piece, reaching
within 1e-6 of every endpoint).

Usage:
  python3 step8_piecewise.py quick    # reruns the failure cases found by step7b / step7c
  python3 step8_piecewise.py sweep    # same sweep as step7_verify quick, piecewise version
"""
import sys
import numpy as np
from scipy.optimize import linprog
from step6_hintz_lp import dual_hintz_form
from step5_window import ctdc_window

P = np.polynomial.polynomial
Leg = np.polynomial.legendre.Legendre


def coeff_arrays(form, zz):
    t, mu, q, L = form['t'], form['mu'], form['q'], form['L']
    return dict(q=P.polyval(zz, q), qp=P.polyval(zz, P.polyder(q)), qpp=P.polyval(zz, P.polyder(q, 2)),
                mu=P.polyval(zz, mu), mup=P.polyval(zz, P.polyder(mu)),
                L=P.polyval(zz, L), Lp=P.polyval(zz, P.polyder(L)))


def CD_rows(cf, f0, f1, f2):
    C = cf['mup']*f0*cf['q'] - 2*cf['mu']*f0*cf['qp'] - cf['mu']*f1*cf['q']
    D = cf['L']*cf['q']*f1 + cf['Lp']*cf['q']*f0 + 0.5*(cf['mup']*(f1*cf['qp'] + f0*cf['qpp'])
                                                       + cf['mu']*(f2*cf['qp'] + 2*f1*cf['qpp']))
    return C, D


def build_basis(a_, b_, K, first_event):
    """Basis polynomials on [a_, b_].  If first_event, each is multiplied by (z - 1)."""
    zm1 = Leg.fit([a_, b_], [a_ - 1.0, b_ - 1.0], 1, domain=[a_, b_])
    out = []
    for k in range(K):
        c = np.zeros(k + 1); c[k] = 1
        p_ = Leg(c, domain=[a_, b_])
        out.append(p_*zm1 if first_event else p_)
    return out


def piecewise_lp(form, z0, side, K=8, npieces=4, N=150, grading=1.0, extra=None):
    t = form['t']
    A, B = (z0, t) if side == 'cosmo' else (1.0, z0)
    # breakpoints, optionally graded towards the horizon end
    s_ = np.linspace(0, 1, npieces + 1)**grading
    if side == 'cosmo':
        bps = B - (B - A)*s_[::-1]  # graded towards t
        bps = np.sort(bps)
    else:
        bps = A + (B - A)*s_        # graded towards 1
    pieces = [(bps[i], bps[i + 1]) for i in range(npieces)]
    nv = K*npieces + 1             # coefficients + tau
    bases = [build_basis(a, b, K, side == 'event' and i == 0) for i, (a, b) in enumerate(pieces)]
    A_ub, b_ub, A_eq, b_eq = [], [], [], []
    # global scales for normalisation (fixed by evaluating C, D magnitudes of a constant-ish f)
    for i, (a, b) in enumerate(pieces):
        zz = np.linspace(a, b, N)
        if i == 0: zz = zz[1:]
        if i == npieces - 1: zz = zz[:-1]
        if extra is not None and i in extra and len(extra[i]):
            zz = np.concatenate([zz, extra[i]])
        cf = {k_: v_[:, None] for k_, v_ in coeff_arrays(form, zz).items()}
        F0 = np.array([p_(zz) for p_ in bases[i]]).T
        F1 = np.array([p_.deriv(1)(zz) for p_ in bases[i]]).T
        F2 = np.array([p_.deriv(2)(zz) for p_ in bases[i]]).T
        RC, RD = CD_rows(cf, F0, F1, F2)
        sC = np.maximum(np.abs(RC).max(axis=1, keepdims=True), 1e-300)
        sD = np.maximum(np.abs(RD).max(axis=1, keepdims=True), 1e-300)
        sF = np.maximum(np.abs(F0).max(axis=1, keepdims=True), 1e-300)
        for R, sc, sign in ((RC, sC, 1.0), (RD, sD, 1.0), (F0, sF, -1.0)):
            blk = np.zeros((R.shape[0], nv))
            blk[:, i*K:(i + 1)*K] = sign*R/sc
            blk[:, -1] = 1.0
            A_ub.extend(blk.tolist()); b_ub.extend([0.0]*R.shape[0])
    # continuity of f and kink condition at interior breakpoints
    for i in range(npieces - 1):
        b = pieces[i][1]
        row = np.zeros(nv)
        row[i*K:(i + 1)*K] = [p_(b) for p_ in bases[i]]
        row[(i + 1)*K:(i + 2)*K] = [-p_(b) for p_ in bases[i + 1]]
        A_eq.append(row.tolist()); b_eq.append(0.0)
        qpb = P.polyval(b, P.polyder(form['q']))
        row = np.zeros(nv)
        row[(i + 1)*K:(i + 2)*K] = [qpb*p_.deriv(1)(b) for p_ in bases[i + 1]]
        row[i*K:(i + 1)*K] = [-qpb*p_.deriv(1)(b) for p_ in bases[i]]
        A_ub.append(row.tolist()); b_ub.append(0.0)    # q'(b)(f'_R - f'_L) <= 0
    # seam condition at z0 and normalisation
    qp0, qpp0 = P.polyval(z0, P.polyder(form['q'])), P.polyval(z0, P.polyder(form['q'], 2))
    if side == 'cosmo':
        i0 = 0; bas = bases[0]
        seam = np.zeros(nv); seam[:K] = [p_.deriv(1)(z0)*qp0 + p_(z0)*qpp0 for p_ in bas]
        A_ub.append(seam.tolist()); b_ub.append(0.0)    # (f q')'(z0) <= 0
        nrm = np.zeros(nv); nrm[(npieces - 1)*K:npieces*K] = [p_(t) for p_ in bases[-1]]
    else:
        bas = bases[-1]
        seam = np.zeros(nv); seam[(npieces - 1)*K:npieces*K] = [-(p_.deriv(1)(z0)*qp0 + p_(z0)*qpp0) for p_ in bas]
        A_ub.append(seam.tolist()); b_ub.append(0.0)    # (f q')'(z0) >= 0
        nrm = np.zeros(nv); nrm[(npieces - 1)*K:npieces*K] = [p_(z0) for p_ in bases[-1]]
    A_eq.append(nrm.tolist()); b_eq.append(1.0)
    bounds = [(-1e4, 1e4)]*(nv - 1) + [(None, 1.0)]
    c = np.zeros(nv); c[-1] = -1.0
    res = linprog(c, A_ub=np.array(A_ub), b_ub=b_ub, A_eq=np.array(A_eq), b_eq=b_eq, bounds=bounds, method='highs')
    if not res.success:
        return None, None
    coef = res.x[:-1].reshape(npieces, K)
    fs = [sum(coef[i][k]*bases[i][k] for k in range(K)) for i in range(npieces)]
    return res.x[-1], (pieces, fs)


def fine_check(form, z0, side, sol, Nf=2000, return_points=False):
    pieces, fs = sol
    worstC = worstD = -np.inf; minf = np.inf; fmax = 0
    bad = {}
    for i, ((a, b), f) in enumerate(zip(pieces, fs)):
        eps = 1e-6*(b - a)
        zz = np.concatenate([np.linspace(a + eps, a + 0.02*(b - a), Nf//4), np.linspace(a + 0.02*(b - a), b - 0.02*(b - a), Nf//2),
                             np.linspace(b - 0.02*(b - a), b - eps, Nf//4)])
        cf = coeff_arrays(form, zz)
        f0, f1, f2 = f(zz), f.deriv(1)(zz), f.deriv(2)(zz)
        C, D = CD_rows(cf, f0, f1, f2)
        # normalise pointwise by the size of the individual terms (relative sign test)
        nC = np.abs(cf['mup']*f0*cf['q']) + np.abs(2*cf['mu']*f0*cf['qp']) + np.abs(cf['mu']*f1*cf['q']) + 1e-300
        nD = (np.abs(cf['L']*cf['q']*f1) + np.abs(cf['Lp']*cf['q']*f0) + 0.5*np.abs(cf['mup']*(f1*cf['qp'] + f0*cf['qpp']))
              + 0.5*np.abs(cf['mu']*(f2*cf['qp'] + 2*f1*cf['qpp'])) + 1e-300)
        worstC = max(worstC, np.max(C/nC)); worstD = max(worstD, np.max(D/nD))
        minf = min(minf, np.min(f0)); fmax = max(fmax, np.max(np.abs(f0)))
        viol = (C/nC > 1e-10) | (D/nD > 1e-10) | (f0 < 0)
        if np.any(viol):
            idx = np.where(viol)[0]
            bad[i] = zz[idx[np.linspace(0, len(idx) - 1, min(len(idx), 40)).astype(int)]]
    out = dict(worstC=worstC, worstD=worstD, minf_rel=minf/fmax)
    return (out, bad) if return_points else out


def run(Lam, a, m, lam, frac, s=-2, M=1.0, configs=((8, 4, 1.0), (8, 8, 2.0), (10, 12, 3.0))):
    lo, hi, _, _ = ctdc_window(Lam, a, M)
    form = dual_hintz_form(Lam, a, M, s, m*(lo + frac*(hi - lo)), m, lam, sign0=-1)
    t = form['t']
    zr = [r.real for r in P.polyroots(form['q']) if abs(r.imag) < 1e-12 and 1 < r.real < t]
    if len(zr) != 1:
        return 'NOZERO', {}
    z0 = zr[0]
    Lq = P.polyval(t, form['L'])*P.polyval(t, form['q'])
    info = dict(z0=z0, t=t, Lq=Lq)
    ok_all = Lq > 0
    for side in ('cosmo', 'event'):
        good = False
        for (K, npc, grad) in configs:
            tau, sol = piecewise_lp(form, z0, side, K, npc, grading=grad)
            if tau is None or tau < 0:
                info[side] = f'lp tau={tau}'
                continue
            fc = fine_check(form, z0, side, sol)
            info[side] = dict(cfg=(K, npc, grad), tau=tau, **fc)
            if fc['worstC'] <= 1e-10 and fc['worstD'] <= 1e-10 and fc['minf_rel'] >= -1e-10:
                good = True; break
        ok_all &= good
    return ('CLOSED' if ok_all else 'OPEN'), info


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    if mode == 'quick':
        cases = [(0.05, 0.6, 1, 2.0, 0.05), (0.05, 0.6, 1, 2.0, 0.95), (0.05, 0.6, 2, 50.0, 0.65),
                 (0.02, 0.9, 1, 2.0, 0.65), (0.02, 0.9, 1, 50.0, 0.95), (0.02, 0.9, 1, 2.0, 0.95),
                 (0.001, 0.95, 2, 10.0, 0.35), (0.001, 0.95, 2, 10.0, 0.65), (0.001, 0.95, 1, 2.0, 0.95)]
    else:
        cases = [(L_, a_, m, lam, fr) for (L_, a_) in [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95), (0.05, 0.9877)]
                 for m in (1, 2, 4) for lam in (2.0, 10.0, 50.0) for fr in (0.05, 0.35, 0.65, 0.95)]
    n = c = 0
    for case in cases:
        st, info = run(*case)
        n += 1; c += st == 'CLOSED'
        def fmt(v):
            if isinstance(v, dict):
                return f"cfg={v['cfg']} tau={v['tau']:.1e} C={v['worstC']:+.1e} D={v['worstD']:+.1e} fmin={v['minf_rel']:+.1e}"
            return str(v)
        if mode == 'quick' or st != 'CLOSED':
            print(f"{st:6} Lam={case[0]} a={case[1]} m={case[2]} lam={case[3]} frac={case[4]}  z0={info.get('z0', 0):.4f} t={info.get('t', 0):.3f} Lq={info.get('Lq', 0):+.2e}")
            print(f"         cosmo: {fmt(info.get('cosmo'))}\n         event: {fmt(info.get('event'))}")
    print(f'\nSUMMARY: {c}/{n} closed')
