"""
Box certification of the event lemma over parameter boxes.

Parameters: t, zf (z0 = 1 + zf (t-1)), A, Bf (B = Bf * Bfloor(t), with Bfloor(t) = 5.5 t - 1.25 proved for lambda >= 1.5).
Variable: tau in [0,1], with z = 1 + tau (z0-1), x = tau X0, s = (1-tau) X0, u = tau, X0 = z0 - 1.
Multiplier: h = g psi + p (1-u)/x, with p = c * min(0.9, pstar), pstar = 2 X0 L(1)/(t-1); (delta, beta, c) fixed per box.
Certify Cc < 0 and x*Dd < 0 for all tau in [0,1] and all parameters in the box, by bisection in tau.
If a box fails, bisect it in its widest relative parameter and retry (up to a depth limit).
"""
from mpmath import iv, mpf
import itertools, time

iv.dps = 25
MULTS = [(0.5, 1.0, 0.5), (1.0, 1.0, 0.5), (0.5, 4.0, 0.5), (0.5, 0.25, 0.5), (2.0, 1.0, 0.2), (0.5, 25.0, 0.2), (1.0, 4.0, 0.05)]


def terms(tau, t, zf, A, Bf, delta, beta, c):
    X0 = zf*(t - 1); z0 = 1 + X0
    B = Bf*(5.5*t - 1.25)
    x = tau*X0; s = (1 - tau)*X0; z = 1 + x; u = tau
    mu_over_x = z*(t - z)
    mup = -3*z**2 + 2*(1 + t)*z - t
    L = A/z + B - 3*z; Lp = -A/z**2 - 3; Lpp = 2*A/z**3
    Lz0 = A/z0 + B - 3*z0; muz0 = z0*X0*(t - z0)
    kappa = Lz0/muz0
    L1 = A + B - 3
    pst = 2*X0*L1/(t - 1)
    p = c*pst
    if p.b > 0.9: p = iv.mpf([min(p.a, 0.9), 0.9]) if p.a < 0.9 else iv.mpf(0.9)
    E = iv.exp(-beta*kappa*s**2)
    psi = (1 + delta)*(1 - E); psip = -(1 + delta)*2*beta*kappa*s*E
    g = -Lp/L; gp = -(Lpp*L - Lp**2)/L**2
    xh = x*g*psi + p*(1 - u)
    x2hh = x**2*(gp*psi + g*psip) - p + xh**2
    H = 2*x*mu_over_x + s*mup
    Cc = s*xh*mu_over_x - H
    xDd = 0.5*(mup*xh + mu_over_x*x2hh) - s*(L*xh + x*Lp)
    # H and the leading terms vanish at tau = 0 or 1 in some factors; work with normalised forms:
    return Cc, xDd, L, p


def cert_box(box, mult, depth_tau=22, max_n=20000):
    t, zf, A, Bf = [iv.mpf(list(b)) for b in box]
    stack = [(mpf(0), mpf(1), 0)]; n = 0
    while stack:
        a, b, d = stack.pop(); n += 1
        if n > max_n: return False
        tau = iv.mpf([a, b])
        Cc, xDd, L, p = terms(tau, t, zf, A, Bf, *mult)
        if L.a > 0 and p.a > 0 and Cc.b < 0 and xDd.b < 0:
            continue
        if d >= depth_tau: return False
        m = (a + b)/2
        stack += [(a, m, d + 1), (m, b, d + 1)]
    return True


def certify_region(box, depth=0, max_depth=10, log=None):
    """Return the list of certified leaf boxes, or the failures."""
    for mult in MULTS:
        if cert_box(box, mult):
            return [(box, mult)], []
    if depth >= max_depth:
        return [], [box]
    # split the widest dimension (relative width; A and Bf in log scale)
    import math
    widths = []
    for i, (lo, hi) in enumerate(box):
        widths.append(math.log(hi/lo) if i in (2, 3) else (hi - lo)/(abs(hi) + abs(lo) + 1e-30))
    i = max(range(4), key=lambda j: widths[j])
    lo, hi = box[i]
    mid = (lo*hi)**0.5 if i in (2, 3) else (lo + hi)/2
    b1 = list(box); b2 = list(box); b1[i] = (lo, mid); b2[i] = (mid, hi)
    ok1, f1 = certify_region(tuple(b1), depth + 1, max_depth)
    ok2, f2 = certify_region(tuple(b2), depth + 1, max_depth)
    return ok1 + ok2, f1 + f2


if __name__ == '__main__':
    import sys
    region = ((2.0, 3.0), (0.2, 0.8), (0.5, 50.0), (1.0, 3.0))
    t0 = time.time()
    leaves, fails = certify_region(region)
    print(f'region {region}: {len(leaves)} certified leaf boxes, {len(fails)} uncertified, {time.time()-t0:.1f}s')
    for f in fails[:10]: print('  uncertified', f)
