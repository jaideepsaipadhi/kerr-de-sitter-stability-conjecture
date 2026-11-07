"""Exact delta_C(sigma) over Z[i] by multimodular CRT with a proved Hadamard bound (item 6, constraint minor).

delta_C = det of the selected 66x66 minor of C(s) = C0 + C1 s + C2 s^2 (rows scaled by their common denominators,
a nonzero constant factor that does not move roots).  Coefficient bound: for |s| = 1,
|det| <= prod_i ||row_i(s)||_2 <= prod_i ( sum_j (sum_k |c_ijk|)^2 )^(1/2) =: B, and every coefficient has |.| <= B
(Cauchy).  For each prime p = 1 mod 4 we get the coefficient images under i -> +iota and i -> -iota, hence x, y mod p
for c = x + i y.  CRT over primes with product > 2B+1 then determines x, y exactly (symmetric residues).
"""
import sys, math, random, time
from math import gcd
from fractions import Fraction as Fr
import flint
from sympy import nextprime
from step55_exceptional import build_point, mod_mat, rank_piv, P as P0, I_MOD

def lcm(a, b): return a*b//gcd(a, b)

def exact_delta(Cq, rows=None, cols=None, verbose=True):
    random.seed(1); s0 = random.randrange(P0)
    r, prow, pcol = rank_piv(mod_mat(Cq, s0, rows, cols))
    rr = [rows[i] for i in prow] if rows is not None else prow
    cc = [cols[j] for j in pcol] if cols is not None else pcol
    # integer-ize rows
    Z = [[[None]*r for _ in range(r)] for _ in range(3)]   # Z[k][a][b] = (xr, xi) ints
    for a, i in enumerate(rr):
        L = 1
        for j in cc:
            for k in range(3):
                for q in Cq[k][i][j]: L = lcm(L, q.denominator)
        for b, j in enumerate(cc):
            for k in range(3):
                x, y = Cq[k][i][j]; Z[k][a][b] = (int(x*L), int(y*L))
    logB = 0.0
    for a in range(r):
        s = sum((sum(math.hypot(*Z[k][a][b]) for k in range(3)))**2 for b in range(r))
        logB += 0.5*math.log2(s) if s > 0 else -1e9
    need = logB + 3
    D = 2*r; xs = list(range(1, D + 2))
    mods, resx, resy = 1, [0]*(D + 1), [0]*(D + 1)
    p = (1 << 62)
    nprimes = 0; t0 = time.time()
    while math.log2(mods) < need:
        p = nextprime(p)
        while p % 4 != 1: p = nextprime(p)
        iota = next(pow(g, (p - 1)//4, p) for g in range(2, 100) if pow(g, (p - 1)//2, p) == p - 1)
        imgs = []
        for io in (iota, p - iota):
            Ms = []
            ys = []
            for x in xs:
                x2 = x*x % p
                M = flint.nmod_mat(r, r, [((Z[0][a][b][0] + io*Z[0][a][b][1]) + (Z[1][a][b][0] + io*Z[1][a][b][1])*x
                                           + (Z[2][a][b][0] + io*Z[2][a][b][1])*x2) % p
                                          for a in range(r) for b in range(r)], p)
                ys.append(int(M.det()))
            V = flint.nmod_mat(D + 1, D + 1, [pow(x, e, p) for x in xs for e in range(D + 1)], p)
            sol = V.solve(flint.nmod_mat(D + 1, 1, ys, p))
            imgs.append([int(sol[e, 0]) for e in range(D + 1)])
        inv2 = pow(2, p - 2, p); inv2i = pow(2*iota % p, p - 2, p)
        xp = [(u + v)*inv2 % p for u, v in zip(*imgs)]; yp = [(u - v)*inv2i % p for u, v in zip(*imgs)]
        # CRT
        newm = mods*p; f = pow(mods, -1, p)
        resx = [(rx + mods*((xv - rx) * f % p)) % newm for rx, xv in zip(resx, xp)]
        resy = [(ry + mods*((yv - ry) * f % p)) % newm for ry, yv in zip(resy, yp)]
        mods = newm; nprimes += 1
    sym = lambda v: v - mods if v > mods//2 else v
    cx = [sym(v) for v in resx]; cy = [sym(v) for v in resy]
    coeffs = list(zip(cx, cy))
    nzk = [k for k, c in enumerate(coeffs) if c != (0, 0)]
    if verbose:
        print(f'    size {r}, log2 B = {logB:.0f}, primes used {nprimes}, {time.time()-t0:.0f}s; '
              f'exact degree {max(nzk)}, exact ord_0 {min(nzk)}', flush=True)
    return coeffs, min(nzk), max(nzk)

if __name__ == '__main__':
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    pts = [(Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)), (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5))]
    for pt in pts:
        B = build_point(*pt, m)
        coeffs, o, d = exact_delta(B['C'])
        lo, hi = coeffs[o], coeffs[d]
        pr = lambda c: (c[0] + I_MOD*c[1]) % P0 != 0   # not in the prime of Z[i] used in step55
        print(f'm={m} pt={pt}: delta_C exact ord {o}, deg {d}; lowest coeff P-unit? {pr(lo)}; leading coeff P-unit? {pr(hi)}', flush=True)
        print(f'    lowest coeff bits {max(abs(lo[0]),abs(lo[1])).bit_length()}, leading coeff bits {max(abs(hi[0]),abs(hi[1])).bit_length()}', flush=True)
