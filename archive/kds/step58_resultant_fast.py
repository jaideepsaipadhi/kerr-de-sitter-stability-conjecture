"""Items 4/5: joint (sigma, m) uniformity via R(mu) = Res_sigma(F1, F2), F_i = E_i / sigma^6, all mod P, fast.

E_i(sigma, mu) (product of the 3 exact minors at test point i) is interpolated ONCE as a bivariate polynomial mod P on
a G x G grid (G = 40 > 29 + slack, top rows/cols checked zero, plus random spot checks).  The FORMAL-degree Sylvester
resultant commutes with reduction mod P, so if R mod P has no root in F_P then for EVERY integer m, R(m) != 0 exactly,
i.e. F1(., m), F2(., m) have no common complex root.
Caveat recorded for the lift (see LEMMAS 68): uses exact divisibility sigma^6 | E_i.
"""
import sys, time, random
import numpy as np
import flint
from fractions import Fraction as Fr
from step56_m_uniform import point_data, eval_mod
from step55_exceptional import P, interp

def fdet(M):
    n = M.shape[0]
    return int(flint.nmod_mat(n, n, [int(x) for x in M.flat], P).det())

def E_eval(minors, s, u):
    v = 1
    for (Mm, rr, cc) in minors: v = v*fdet(eval_mod(Mm, s, u, rows=rr, cols=cc)) % P
    return v

def biv(minors, G):
    xs = list(range(1, G + 1))
    rows = [interp(xs, [E_eval(minors, s, u) for u in xs]) for s in xs]          # in mu, for each sigma
    rows = [(r + [0]*G)[:G] for r in rows]
    cols = [interp(xs, [rows[i][l] for i in range(G)]) for l in range(G)]       # in sigma, for each mu-power l
    C = [[(cols[l] + [0]*G)[k] for l in range(G)] for k in range(G)]             # C[k][l]: sigma^k mu^l
    return C

def ev2(C, s, u):
    return sum(C[k][l]*pow(s, k, P) % P*pow(u, l, P) for k in range(len(C)) for l in range(len(C)) if C[k][l]) % P

def degs(C):
    ks = [k for k in range(len(C)) if any(C[k])]; ls = [l for l in range(len(C)) if any(C[k][l] for k in range(len(C)))]
    return min(ks), max(ks), max(ls)

if __name__ == '__main__':
    t0 = time.time(); G = 40
    p1 = (Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)); p2 = (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5))
    Cs = []
    for p in (p1, p2):
        M, I = point_data(*p, t0)
        C = biv(M, G); o, dk, dl = degs(C)
        random.seed(11); ok = all(ev2(C, s, u) == E_eval(M, s, u) for s, u in [(random.randrange(P), random.randrange(P)) for _ in range(3)])
        print(f'  E at {p}: sigma-order {o}, deg_sigma {dk}, deg_mu {dl} (grid {G}); random spot checks OK: {ok}  [{time.time()-t0:.0f}s]', flush=True)
        Cs.append((C, o, dk, dl))
    (C1, o1, d1, l1), (C2, o2, d2, l2) = Cs
    F1 = [C1[k] for k in range(o1, d1 + 1)]; F2 = [C2[k] for k in range(o2, d2 + 1)]
    n1, n2 = len(F1) - 1, len(F2) - 1                         # formal sigma-degrees
    DR = n2*l1 + n1*l2
    def sylv(u):
        a = [sum(c*pow(u, l, P) for l, c in enumerate(row)) % P for row in F1]
        b = [sum(c*pow(u, l, P) for l, c in enumerate(row)) % P for row in F2]
        n = n1 + n2; M = np.zeros((n, n), dtype=np.int64)
        for i in range(n2):
            for j, c in enumerate(reversed(a)): M[i, i + j] = c
        for i in range(n1):
            for j, c in enumerate(reversed(b)): M[n2 + i, i + j] = c
        return fdet(M)
    us = list(range(1, DR + 2)); R = interp(us, [sylv(u) for u in us])
    Rp = flint.nmod_poly(R, P)
    chk = all(int(Rp(x)) == sylv(x) for x in (DR + 5, DR + 77))
    print(f'R(mu): degree {Rp.degree()} (bound {DR}), interpolant check {chk}, identically zero {Rp.is_zero()}  [{time.time()-t0:.0f}s]', flush=True)
    # roots in F_P:  gcd(R, mu^P - mu)
    x = flint.nmod_poly([0, 1], P); acc = flint.nmod_poly([1], P); base = x; e = P
    while e:
        if e & 1: acc = (acc*base) % Rp
        base = (base*base) % Rp; e >>= 1
    g = Rp.gcd(acc - x)
    print(f'deg gcd(R, mu^P - mu) = {g.degree()}  => R has a root in F_P: {g.degree() > 0}', flush=True)
    if g.degree() > 0:
        rts = [r for r in range(-2000, 2001) if int(Rp(r % P)) == 0]
        print('  small integer mu with R(mu) = 0 mod P:', rts)
        print('  F_P roots:', g.factor())
