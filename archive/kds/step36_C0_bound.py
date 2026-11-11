"""Step 36: prove C0 > (2t-1)/(3(1-alpha)).
C0 = 3/(Lam (r1-r0)(r2-r3)),  t = (r1-r3)(r2-r0)/((r1-r0)(r2-r3)),  3(1-alpha)/Lam = -e2.
Inequality <=> F = -3 e2 - 2(r1-r3)(r2-r0) + (r1-r0)(r2-r3) > 0, with r3 = -(r0+r1+r2), 0 < r0 < r1 < r2."""
import sympy as sp
x, y, w = sp.symbols('r0 r1 r2', positive=True)
r3 = -(x + y + w)
rs = [x, y, w, r3]
e2 = sum(rs[i]*rs[j] for i in range(4) for j in range(i + 1, 4))
F = sp.expand(-3*e2 - 2*(y - r3)*(w - x) + (y - x)*(w - r3))
print('F =', sp.factor(F))
# substitute r1 = r0 + u, r2 = r1 + v (u, v > 0) and look at the signs of the monomials
u, v = sp.symbols('u v', positive=True)
G = sp.expand(F.subs({y: x + u, w: x + u + v}))
print('F(r0, r0+u, r0+u+v) =', G)
print('all coefficients positive:', all(c > 0 for c in sp.Poly(G, x, u, v).coeffs()))
# numerical check of C0 against the formula
import numpy as np
from step3_dual import kds_roots
for (Lam, a) in [(0.05, 0.6), (0.001, 0.95), (0.01, 0.99)]:
    r3n, r0n, r1n, r2n = [float(q) for q in kds_roots(Lam, a, 1.0)]
    print(f'Lam={Lam} a={a}: C0 formula = {3/(Lam*(r1n - r0n)*(r2n - r3n)):.6f}')
