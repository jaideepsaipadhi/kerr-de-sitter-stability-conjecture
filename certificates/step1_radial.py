"""
Step 1: the spin-s radial Teukolsky operator on Kerr-de Sitter, and checks on its conventions.

Radial equation (Suzuki-Takasugi-Umetsu, gr-qc/9905040, eq. (2.8) with Q=0):

  Delta^{-s} (Delta^{s+1} R')' + [ Xi^2 K^2 - i s Xi K Delta' ] / Delta
      + 4 i s Xi omega r - (2 Lambda/3) c2 r^2 + s(1 - alpha) - lam = 0

  Delta = (r^2+a^2)(1 - Lambda r^2/3) - 2 M r,  K = omega (r^2+a^2) - a m,
  alpha = Lambda a^2/3,  Xi = 1 + alpha.

The paper prints two versions of the r^2 coefficient: c2 = (s+1)(2s+1) and
c2 = (s-1)(2s-1).  Check (A) settles which one is right: only the correct
coefficient makes r = infinity a removable singularity.  Infinity is not a
point of the geometry, so it has to be removable.

Checks:
  (A) Indicial exponents at r = infinity for both versions of c2.
  (B) Lambda -> 0 gives the Kerr Teukolsky radial operator.
"""
import sympy as sp

r, a, M, Lam, w, s, lam, m = sp.symbols('r a M Lambda omega s lambda m')
rho = sp.symbols('rho')


def radial_coeffs(c2):
    alpha = Lam*a**2/3
    Xi = 1 + alpha
    D = (r**2 + a**2)*(1 - Lam*r**2/3) - 2*M*r
    K = w*(r**2 + a**2) - a*m
    # Multiply the equation by Delta: Delta R'' + (s+1) Delta' R' + Q R = 0
    Q = (Xi**2*K**2 - sp.I*s*Xi*K*sp.diff(D, r))/D + 4*sp.I*s*Xi*w*r \
        - sp.Rational(2, 3)*Lam*c2*r**2 + s*(1 - alpha) - lam
    p = (s + 1)*sp.diff(D, r)/D
    q = Q/D
    return D, p, q


def indicial_at_infinity(c2):
    D, p, q = radial_coeffs(c2)
    # Look for R ~ r^rho as r -> infinity.
    # Leading behaviour as r -> infinity via leading coefficients of rational
    # functions (sympy.limit is far too slow on these expressions).
    def lead(expr):
        num, den = sp.fraction(sp.together(sp.expand(expr)))
        num, den = sp.Poly(sp.expand(num), r), sp.Poly(sp.expand(den), r)
        if num.degree() < den.degree():
            return sp.Integer(0)
        assert num.degree() == den.degree(), 'unexpected growth at infinity'
        return sp.simplify(num.LC()/den.LC())
    p0 = lead(p*r)
    q0 = lead(q*r**2)
    return sp.factor(sp.expand(rho*(rho - 1) + p0*rho + q0)), p0, q0


if __name__ == '__main__':
    for name, c2 in [('(s+1)(2s+1)', (s + 1)*(2*s + 1)), ('(s-1)(2s-1)', (s - 1)*(2*s - 1))]:
        ind, p0, q0 = indicial_at_infinity(c2)
        print(f'c2 = {name}: indicial polynomial at r = inf is {ind}; roots {sp.solve(ind, rho)}')

    # (B) Kerr limit.  Teukolsky's separation constant lam_T corresponds to lam = lam_T + s.
    D, p, q = radial_coeffs((s + 1)*(2*s + 1))
    Qkerr = sp.simplify(sp.together(q*D).subs(Lam, 0))
    Dk = r**2 - 2*M*r + a**2
    Kk = w*(r**2 + a**2) - a*m
    lamT = sp.symbols('lambda_T')
    teuk = (Kk**2 - 2*sp.I*s*(r - M)*Kk)/Dk + 4*sp.I*s*w*r - lamT
    print('Kerr limit matches Teukolsky:', sp.simplify(Qkerr.subs(lam, lamT + s) - teuk) == 0)
