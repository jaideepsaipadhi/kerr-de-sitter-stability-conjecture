import time, flint
from fractions import Fraction as Fr
from step56_m_uniform import point_data
from step55_exceptional import P
from step58_resultant_fast import biv, degs
t0 = time.time(); G = 81
pts = [(Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)), (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5)), (Fr(1, 20), Fr(3, 5), Fr(4), Fr(1, 3))]
Fs = []
for p in pts:
    M, I = point_data(*p, t0); C = biv(M, G); o, dk, dl = degs(C)
    F = [C[k] for k in range(o, dk + 1)]            # F[k][l]: sigma^k mu^l
    Fs.append(F)
    ordmu = min(l for k in range(len(F)) for l in range(len(F[0])) if F[k][l])
    f0 = flint.nmod_poly([F[k][0] for k in range(len(F))], P)
    print(f'pt {p}: mu-order of F = {ordmu}; F(sigma,0) has sigma-degree {f0.degree()}, factor degrees {[ (f.degree(), e) for f, e in f0.factor()[1]]}', flush=True)
f = [flint.nmod_poly([F[k][0] for k in range(len(F))], P) for F in Fs]
print('gcd F1(.,0),F2(.,0):', f[0].gcd(f[1]).degree(), ' gcd all three:', f[0].gcd(f[1]).gcd(f[2]).degree())
# sigma-only structure: does F depend on mu only via mu^2 ?  and via sigma*mu combos?
for i, F in enumerate(Fs):
    odd = any(F[k][l] for k in range(len(F)) for l in range(1, len(F[0]), 2))
    print(f'pt{i}: odd powers of mu present: {odd}')
import numpy as np
from step58_resultant_fast import fdet
from step55_exceptional import interp
Hs = [[row[6:] for row in F] for F in Fs]
h0 = [flint.nmod_poly([H[k][0] for k in range(len(H))], P) for H in Hs]
for i, h in enumerate(h0): print(f'H{i}(sigma,0): degree {h.degree()}, factors {[(f.degree(), e) for f, e in h.factor()[1]]}')
print('gcd H1(.,0),H2(.,0):', h0[0].gcd(h0[1]).degree(), ' all three:', h0[0].gcd(h0[1]).gcd(h0[2]).degree())
def R(i, j):
    F1, F2 = Hs[i], Hs[j]; n1, n2 = len(F1) - 1, len(F2) - 1
    l1 = max(l for r in F1 for l, c in enumerate(r) if c); l2 = max(l for r in F2 for l, c in enumerate(r) if c); DR = n2*l1 + n1*l2
    def sylv(u):
        a = [sum(c*pow(u, l, P) for l, c in enumerate(r)) % P for r in F1]; b = [sum(c*pow(u, l, P) for l, c in enumerate(r)) % P for r in F2]
        n = n1 + n2; Mx = np.zeros((n, n), dtype=np.int64)
        for k in range(n2):
            for jj, c in enumerate(reversed(a)): Mx[k, k + jj] = c
        for k in range(n1):
            for jj, c in enumerate(reversed(b)): Mx[n2 + k, k + jj] = c
        return fdet(Mx)
    us = list(range(1, DR + 2)); return flint.nmod_poly(interp(us, [sylv(u) for u in us]), P)
Rs = [R(0, 1), R(0, 2), R(1, 2)]
for f in Rs:
    c = f.coeffs(); print('R(H) degree', f.degree(), 'order at 0', next(k for k, x in enumerate(c) if int(x)))
g = Rs[0].gcd(Rs[1]).gcd(Rs[2]); print('gcd of the three R(H):', g)
from math import isqrt
from fractions import Fraction
g9 = h0[0].gcd(h0[1]).gcd(h0[2])
r = [f for f, e in g9.factor()[1]][0]; c = r.coeffs(); root = (-int(c[0])*pow(int(c[1]), P - 2, P)) % P
from step55_exceptional import I_MOD
def ratrec(x):
    B = isqrt(P//2); r0, r1, s0, s1 = P, x % P, 0, 1
    while r1 > B: q = r0//r1; r0, r1 = r1, r0 - q*r1; s0, s1 = s1, s0 - q*s1
    return Fraction(r1, s1) if s1 and abs(s1) <= B else None
print('sigma* mod P =', root, ' as rational:', ratrec(root), ' sigma*/i as rational:', ratrec(root*pow(I_MOD, P - 2, P) % P))
