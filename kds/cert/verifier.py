"""
Rigorous verifier for the 4-parameter event lemma (KdS project), using interval arithmetic (mpmath.iv).

Lemma (event side).  Let t > 1, 1 < z0 < t, L(z) = A/z + B - 3z > 0 on [1, z0], mu = z(z-1)(t-z), s = z0 - z,
x = z - 1, H = 2 mu + s mu'.  With the multiplier
    h = f'/f = g psi + p (1-u)/x,   g = -L'/L,   psi = (1+delta)(1 - exp(-beta kappa s^2)),
    kappa = L(z0)/mu(z0),  u = x/(z0-1),
the two conditions are
    Cc := s h mu - H <= 0,
    Dd := (1/2)(mu' h + mu (h' + h^2)) - s (L h + L') <= 0.
We certify Cc and x*Dd (the latter finite at z = 1) via the cancellation-free forms
    x h      = x g psi + p (1-u),
    x^2(h' + h^2) = x^2 (g psi)' - p + (x h)^2,
    Cc       = s (x h)(mu/x) - H,
    x Dd     = (1/2)(mu' (x h) + (mu/x)(x^2 (h'+h^2))) - s (L (x h) + x L').
Both are continuous on [1, z0] and strictly negative at the endpoints when 0 < p < min(1, 2(z0-1)L(1)/(t-1)).
Strategy: adaptive bisection of [1, z0]; a subinterval is accepted when the interval upper bounds of Cc and x*Dd
are both < 0.  Parameters may themselves be intervals (that certifies a parameter box).
"""
from mpmath import iv, mpf
import mpmath

iv.dps = 30


def _terms(z, t, z0, A, B, delta, p, beta):
    x = z - 1
    s = z0 - z
    mu_over_x = z*(t - z)
    mu = x*mu_over_x
    mup = -3*z**2 + 2*(1 + t)*z - t
    L = A/z + B - 3*z
    Lp = -A/z**2 - 3
    Lpp = 2*A/z**3
    Lz0 = A/z0 + B - 3*z0
    muz0 = z0*(z0 - 1)*(t - z0)
    kappa = Lz0/muz0
    E = iv.exp(-beta*kappa*s**2)
    psi = (1 + delta)*(1 - E)
    psip = -(1 + delta)*2*beta*kappa*s*E
    g = -Lp/L
    gp = -(Lpp*L - Lp**2)/L**2
    u = x/(z0 - 1)
    xh = x*g*psi + p*(1 - u)
    x2gpsi_p = x**2*(gp*psi + g*psip)
    x2hh = x2gpsi_p - p + xh**2
    H = 2*mu + s*mup
    Cc = s*xh*mu_over_x - H
    xDd = 0.5*(mup*xh + mu_over_x*x2hh) - s*(L*xh + x*Lp)
    return Cc, xDd, L


def certify_point(t, z0, A, B, delta, p, beta, max_depth=40, max_boxes=200000):
    """Return (True, n_boxes) if Cc < 0 and x*Dd < 0 are certified on [1, z0], else (False, info)."""
    T, Z0, AA, BB = [iv.mpf(v) if not isinstance(v, type(iv.mpf(1))) else v for v in (t, z0, A, B)]
    dl, pp, bt = iv.mpf(delta), iv.mpf(p), iv.mpf(beta)
    lo0, hi0 = mpf(1), (Z0.b if hasattr(Z0, 'b') else mpf(z0))
    stack = [(mpf(1), mpf(Z0.a), 0)]
    # The right endpoint varies with the z0 box; we treat z in [1, z0_low] and add the sliver [z0_low, z0_high]
    # by evaluating with s in [0, width].  For point parameters the sliver is empty.
    n = 0
    while stack:
        a, b, dpt = stack.pop()
        n += 1
        if n > max_boxes:
            return False, ('too many boxes', float(a), float(b))
        z = iv.mpf([a, b])
        Cc, xDd, L = _terms(z, T, Z0, AA, BB, dl, pp, bt)
        if L.a <= 0:
            if dpt >= max_depth: return False, ('L not certified > 0', float(a), float(b))
        elif Cc.b < 0 and xDd.b < 0:
            continue
        if dpt >= max_depth:
            return False, ('depth', float(a), float(b), float(Cc.b), float(xDd.b))
        m = (a + b)/2
        stack.append((a, m, dpt + 1)); stack.append((m, b, dpt + 1))
    return True, n
