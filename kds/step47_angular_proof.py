"""
Step 47: rigorous angular bound lam >= 3/2 for |m| = 1, 2 (|m| >= 3 is analytic, see LEMMAS 45).
H >= H0 + E (form sense), E = N/(1+al x^2) - s al + 2 al x^2 with
  N = -s^2 al (1-x^2) + (1+al)^2 (xi^2 (1-x^2) + 2 xi (m + s x)) - 2 s al m (1+al) x.
Since 1/(1+al x^2) = 1 - al x^2/(1+al x^2), E = P + R with P = N - s al + 2 al x^2 (polynomial in x) and
|R| <= al sup|N| <= al * C, where C = 4 al + (1+al)^2 (xi^2 + 2 xi (|m|+2)) + 4 al |m| (1+al).
Two-level bound (H0 levels 2 and 8, psi0 = l=2 state, |psi0|^2 ∝ (1-x)^{|m+s|}(1+x)^{|m-s|}):
  lam_min(H0 + P) >= mineig [[2 + <P>, sigma], [sigma, 8 + Pmin]],  sigma^2 = <P^2> - <P>^2,
with Pmin a rigorous lower bound for inf_x P.  So  lam >= that - al C.
<P>, <P^2> are exact polynomials in (al, xi) (Beta moments).  Check >= 3/2 on al in [0, 0.0718], xi in [0, 0.45|m|]
(WLOG m > 0, xi >= 0 by the symmetry (x, m, xi) -> (-x, -m, -xi)), by interval arithmetic on a box grid.
"""
import sympy as sp
from mpmath import iv, mpf

x, al, xi = sp.symbols('x al xi')

def moments_poly(m, s):
    N = -s**2*al*(1 - x**2) + (1 + al)**2*(xi**2*(1 - x**2) + 2*xi*(m + s*x)) - 2*s*al*m*(1 + al)*x
    P = sp.expand(N - s*al + 2*al*x**2)
    w = (1 - x)**abs(m + s)*(1 + x)**abs(m - s)
    Z = sp.integrate(w, (x, -1, 1))
    mom = lambda f: sp.expand(sp.integrate(sp.expand(w*f), (x, -1, 1))/Z)
    p1 = mom(P); p2 = mom(P**2)
    C = 4*al + (1 + al)**2*(xi**2 + 2*xi*(abs(m) + 2)) + 4*al*abs(m)*(1 + al)
    # rigorous lower bound of P on [-1,1]: P = c0 + c1 x + c2 x^2 (degree 2 in x); bound termwise
    c = sp.Poly(P, x).all_coeffs()[::-1] + [0, 0, 0]
    Pmin = c[0] - sp.Abs(c[1]) - sp.Abs(c[2])          # >= inf over |x|<=1 (crude, valid)
    return p1, p2, C, Pmin

def check(m, s, nA=24, nX=120, amax=0.0718, xmax_fac=1.0):
    p1, p2, C, Pmin = moments_poly(m, s)
    f1, f2, fC, fP = [sp.lambdify((al, xi), e, modules=[{'Abs': abs, 'mpf': iv.mpf}, 'mpmath']) for e in (p1, p2, C, Pmin)]
    worst = None
    XM = xmax_fac*abs(m)
    for i in range(nA):
        A = iv.mpf([amax*i/nA, amax*(i + 1)/nA])
        for j in range(nX):
            X = iv.mpf([XM*j/nX, XM*(j + 1)/nX])
            P1, P2, CC, PM = f1(A, X), f2(A, X), fC(A, X), fP(A, X)
            var = P2 - P1*P1
            var = iv.mpf([max(mpf(0), var.a), var.b])
            a11 = 2 + P1; a22 = 8 + PM
            # smaller eigenvalue of [[a11, sig],[sig, a22]] >= (a11 + a22)/2 - sqrt(((a11-a22)/2)^2 + var); evaluate lower end
            lam = (a11 + a22)/2 - iv.sqrt(((a11 - a22)/2)**2 + var) - A*CC
            lo = lam.a
            if worst is None or lo < worst[0]: worst = (lo, float(A.a), float(X.a))
    return worst

if __name__ == '__main__':
    ok = True
    for m in (1, 2):
        for s in (-2, 2):
            w = check(m, s)
            print(f'm={m} s={s:+d}: certified lower bound on lam over the box: {float(w[0]):.4f} (at al~{w[1]:.4f}, xi~{w[2]:.3f})', flush=True)
            ok &= w[0] >= 1.5
    print('ALL >= 3/2:', ok)
