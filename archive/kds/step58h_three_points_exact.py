"""Route (A): gcd(R12, R13) mod p with a third test point.  Any integer m with a sigma exceptional at all three points is
a root of both resultants, hence (mod p) of g_p = gcd(R12, R13)."""
import time, random, numpy as np, flint
from fractions import Fraction as Fr
from step56_m_uniform import point_data, eval_mod
from step55_exceptional import P, interp
from step58_resultant_fast import biv, degs, fdet
t0 = time.time(); G = 81   # exact: deg_sigma, deg_mu <= 80 proved (step58g)
pts = [(Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)), (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5)), (Fr(1, 20), Fr(3, 5), Fr(4), Fr(1, 3))]
Fs = []
for p in pts:
    M, I = point_data(*p, t0); C = biv(M, G); o, dk, dl = degs(C)
    print(f'  {p}: sigma-order {o}, deg_sigma {dk}, deg_mu {dl}', flush=True)
    Fs.append(([C[k] for k in range(o, dk + 1)], dl))
def R(i, j):
    (F1, l1), (F2, l2) = Fs[i], Fs[j]; n1, n2 = len(F1) - 1, len(F2) - 1; DR = n2*l1 + n1*l2
    def sylv(u):
        a = [sum(c*pow(u, l, P) for l, c in enumerate(r)) % P for r in F1]; b = [sum(c*pow(u, l, P) for l, c in enumerate(r)) % P for r in F2]
        n = n1 + n2; Mx = np.zeros((n, n), dtype=np.int64)
        for k in range(n2):
            for jj, c in enumerate(reversed(a)): Mx[k, k + jj] = c
        for k in range(n1):
            for jj, c in enumerate(reversed(b)): Mx[n2 + k, k + jj] = c
        return fdet(Mx)
    us = list(range(1, DR + 2)); return flint.nmod_poly(interp(us, [sylv(u) for u in us]), P)
R12, R13, R23 = R(0, 1), R(0, 2), R(1, 2)
x = flint.nmod_poly([0, 1], P)
def order0(f):
    k = 0
    while int(f.coeffs()[k]) == 0: k += 1
    return k
for nm, f in (('R12', R12), ('R13', R13), ('R23', R23)): print(f'{nm}: degree {f.degree()}, order at mu=0: {order0(f)}', flush=True)
g = R12.gcd(R13).gcd(R23)
print(f'gcd(R12, R13, R23) mod p = {g}   (degree {g.degree()})  [{time.time()-t0:.0f}s]')
