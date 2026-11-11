"""
Box certification, version 2: center-fixed multiplier, with (A,B) handled by the 4 corners of the box.

Within a box, f is fixed as a function of tau = (z-1)/X0, built from the box-center parameters (tc, zfc, Ac, Bc).
Then  x h = x_c h_c(z_c(tau))  and  x^2 h' = x_c^2 h_c'(z_c(tau))  are center quantities.  This is exact: it is the
chain rule for f(z) = F(tau(z)).
   Cc   = s (xh)(mu/x) - H                             : independent of (A,B)
   x Dd = (1/2)(mu'(xh) + (mu/x)(x^2 h' + (xh)^2)) - s (L (xh) + x L'),   L = A/z + B - 3z : LINEAR in (A,B)
For each tau-interval and the actual (t, zf) intervals, x*Dd < 0 at the 4 (A,B) corners gives x*Dd < 0 on the
whole (A,B) rectangle.  L > 0 likewise follows from the corners (L is linear in (A,B)), together with L > 0 at the
smallest A,B (L is decreasing in z, so it suffices to check z = z0).
Endpoint signs: at tau = 0 (z = 1), Cc = X0 (t-1)(p_c X0/X0c ...) -- handled automatically by the interval
evaluation of the cancellation-free forms.
"""
from mpmath import iv, mpf
import math, time

iv.dps = 25
MULTS = [(0.5, 1.0, 0.5), (0.1, 0.05, 0.5), (1.0, 1.0, 0.5), (0.25, 0.25, 0.5), (0.5, 4.0, 0.5), (0.1, 0.25, 0.5), (0.5, 0.25, 0.5),
         (2.0, 1.0, 0.2), (0.5, 25.0, 0.2), (1.0, 4.0, 0.05), (0.1, 0.05, 0.2)]


def Bfloor(t):
    return 5.5*t - 1.25


def center_funcs(tau, tc, zfc, Ac, Bc, delta, beta, c):
    """Return (x h, x^2 h' + (x h)^2) as functions of the tau interval, built from the center parameters."""
    X0 = zfc*(tc - 1); z0 = 1 + X0
    x = tau*X0; s = (1 - tau)*X0; z = 1 + x
    L = Ac/z + Bc - 3*z; Lp = -Ac/z**2 - 3; Lpp = 2*Ac/z**3
    Lz0 = Ac/z0 + Bc - 3*z0; muz0 = z0*X0*(tc - z0)
    kappa = Lz0/muz0
    L1 = Ac + Bc - 3
    p = min(0.9, c*2*X0*L1/(tc - 1))
    E = iv.exp(-beta*kappa*s**2)
    psi = (1 + delta)*(1 - E); psip = -(1 + delta)*2*beta*kappa*s*E
    g = -Lp/L; gp = -(Lpp*L - Lp**2)/L**2
    xh = x*g*psi + p*(1 - tau)
    x2hp = x**2*(gp*psi + g*psip) - p
    return xh, x2hp + xh**2, p


def cert_box(box, mult, depth_tau=24, max_n=40000):
    (t0, t1), (z0a, z0b), (A0, A1), (Bf0, Bf1) = box
    tc, zfc = (t0 + t1)/2, (z0a + z0b)/2
    Ac, Bfc = math.sqrt(A0*A1), math.sqrt(Bf0*Bf1)
    Bc = Bfc*Bfloor(tc)
    t = iv.mpf([t0, t1]); zf = iv.mpf([z0a, z0b])
    corners = [(A, Bf) for A in (A0, A1) for Bf in (Bf0, Bf1)]
    stack = [(mpf(0), mpf(1), 0)]; n = 0
    while stack:
        a, b, d = stack.pop(); n += 1
        if n > max_n: return False
        tau = iv.mpf([a, b])
        xh, X2, p = center_funcs(tau, tc, zfc, Ac, Bc, *mult)
        X0 = zf*(t - 1); x = tau*X0; s = (1 - tau)*X0; z = 1 + x
        mu_over_x = z*(t - 1)*(1 - tau*zf)   # t - z = (t-1)(1 - tau zf), no cancellation
        mup = mu_over_x - x*(2*z - t)
        # C = s[(xh - 1) mu/x + x(2z - t)] - 2 x mu/x, using mu' = mu/x - x(2z - t)  (dependency-safe form)
        Cc = s*((xh - 1)*mu_over_x + x*(2*z - t)) - 2*x*mu_over_x
        ok = Cc.b < 0
        if ok:
            for (A, Bf) in corners:
                B = Bf*Bfloor(t)
                L = A/z + B - 3*z; Lp = -A/z**2 - 3
                xDd = 0.5*(mu_over_x*(xh + X2) - x*(2*z - t)*xh) - s*(L*xh + x*Lp)
                if not (xDd.b < 0 and L.a > 0):
                    ok = False; break
        if ok:
            continue
        if d >= depth_tau: return False
        m = (a + b)/2
        stack += [(a, m, d + 1), (m, b, d + 1)]
    return True


def _fwd(i, v):
    import math
    if i == 0: return math.log(v - 1)                 # t     -> log(t-1)
    if i == 1: return math.log(v/(1 - v))             # zf    -> logit
    return math.log(v)                                # A, Bf -> log


def _inv(i, y):
    import math
    if i == 0: return 1 + math.exp(y)
    if i == 1: return 1/(1 + math.exp(-y))
    return math.exp(y)


def certify_region(region, max_depth=18, max_boxes=4000):
    """Adaptive splitting in corner-adapted coordinates (log(t-1), logit(zf), log A, log Bf).
    Split the dimension that is widest relative to the region's own extent.  Returns (n_certified, failures)."""
    base = [max(_fwd(i, hi) - _fwd(i, lo), 1e-12) for i, (lo, hi) in enumerate(region)]
    todo = [(region, 0)]; done = 0; fails = []; tried = 0
    while todo:
        box, d = todo.pop(); tried += 1
        if tried > max_boxes:
            fails.append(('BOX CAP', box)); fails += [bx for bx, _ in todo]; break
        if any(cert_box(box, m) for m in MULTS):
            done += 1; continue
        if d >= max_depth:
            fails.append(box); continue
        w = [(_fwd(i, hi) - _fwd(i, lo))/base[i] for i, (lo, hi) in enumerate(box)]
        i = max(range(4), key=lambda j: w[j])
        lo, hi = box[i]
        mid = _inv(i, (_fwd(i, lo) + _fwd(i, hi))/2)
        b1 = list(box); b2 = list(box); b1[i] = (lo, mid); b2[i] = (mid, hi)
        todo += [(tuple(b1), d + 1), (tuple(b2), d + 1)]
    return done, fails


if __name__ == '__main__':
    import sys
    region = ((2.0, 3.0), (0.2, 0.8), (0.5, 50.0), (1.0, 3.0))
    t0 = time.time()
    done, fails = certify_region(region)
    print(f'region {region}: {done} boxes certified, {len(fails)} uncertified, {time.time()-t0:.1f}s')
    for f in fails[:5]: print('  uncertified', f)
