"""Sharper C0 bound: 2(1-alpha) C0 >= 2t - 1  <=>  G = -2 e2 - 2(r1-r3)(r2-r0) + (r1-r0)(r2-r3) >= 0."""
import sympy as sp
x, y, w, u, v = sp.symbols('r0 r1 r2 u v', positive=True)
r3 = -(x + y + w); rs = [x, y, w, r3]
e2 = sum(rs[i]*rs[j] for i in range(4) for j in range(i + 1, 4))
G = sp.expand(-2*e2 - 2*(y - r3)*(w - x) + (y - x)*(w - r3))
print('G =', sp.factor(G))
Gs = sp.expand(G.subs({y: x + u, w: x + u + v}))
print('G(r0, r0+u, r0+u+v) =', Gs, '  all coefficients >= 0:', all(c >= 0 for c in sp.Poly(Gs, x, u, v).coeffs()))
