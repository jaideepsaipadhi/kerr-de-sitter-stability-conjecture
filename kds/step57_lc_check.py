"""Item 6: lift the mod-P exceptional-sigma certificate (step55) to an exact statement.

Lemma (Gauss + reduction).  Let E1, E2 in Z[i][s] (denominators cleared), P a prime of Z[i].
Suppose (a) deg E1 = deg(E1 mod P)            [leading coeff not divisible by P]
        (b) ord_0 E1 = ord_0(E1 mod P) =: b    [lowest nonzero coefficient not divisible by P]
        (c) gcd(E1 mod P, E2 mod P) = u s^c.
Then every common complex root of E1, E2 is s = 0.
Proof.  Let h be a primitive gcd of E1, E2 in Z[i][s].  lc(h) | lc(E1), so by (a) h has all its roots P-integral
and h mod P | u s^c, so every root s0 of h satisfies s0 = 0 mod P.  Write E1 = s^b F, F(0) != 0; by (b) F(0) is a
P-unit.  If s0 != 0 is a root of h, then F(s0) = 0, but F(s0) = F(0) mod P (s0 is P-integral and = 0 mod P),
a unit: contradiction.  QED
E1 is a product of minors, so (a),(b) hold for E1 iff they hold for every factor minor.

Exact certificate of (a),(b) per minor: with d_ij = deg entry, v_ij = ord_0 entry, the true degree is
<= max-weight assignment of d over the minor (the tropical bound), and the true order is >= min-weight assignment
of v.  If max-assignment(d) = deg(minor mod P) and min-assignment(v) = ord(minor mod P), then the true degree and
order are pinned (true deg >= mod-P deg, true ord <= mod-P ord, always), i.e. (a),(b) hold exactly.
"""
import sys, time, random
import numpy as np
from scipy.optimize import linear_sum_assignment
from step55_exceptional import build_point, mod_mat, rank_piv, det_mod, interp, P, Fr

def nz(c): return c[0] != 0 or c[1] != 0

def minor_data(Cq, rows=None, cols=None):
    random.seed(1); s0 = random.randrange(P)
    r, prow, pcol = rank_piv(mod_mat(Cq, s0, rows, cols))
    rr = [rows[i] for i in prow] if rows is not None else prow
    cc = [cols[j] for j in pcol] if cols is not None else pcol
    D = 2*r; xs = list(range(1, D + 2))
    poly = interp(xs, [det_mod(mod_mat(Cq, x, rr, cc)) for x in xs])
    degP = len(poly) - 1; ordP = next(k for k, c in enumerate(poly) if c)
    NEG = -10**6
    dm = np.full((r, r), NEG, dtype=float); vm = np.full((r, r), 10**6, dtype=float)
    for a, i in enumerate(rr):
        for b, j in enumerate(cc):
            ks = [k for k in range(3) if nz(Cq[k][i][j])]
            if ks: dm[a, b] = max(ks); vm[a, b] = min(ks)
    ri, ci = linear_sum_assignment(dm, maximize=True); dmax = dm[ri, ci].sum()
    ri, ci = linear_sum_assignment(vm); vmin = vm[ri, ci].sum()
    return r, degP, ordP, int(dmax), int(vmin)

if __name__ == '__main__':
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    pts = [(Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)), (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5))]
    allok = True
    for pt in pts:
        t0 = time.time(); B = build_point(*pt, m)
        mins = [('C', minor_data(B['C']))] + [(f'G{J}', minor_data(B['G'], cols=B['low'][J])) for J in (0, 1)]
        for name, (r, dP, oP, dmax, vmin) in mins:
            ok = (dmax == dP) and (vmin == oP); allok &= ok
            print(f'm={m} pt={pt} {name}: size {r}  deg mod P {dP} vs tropical max {dmax};  '
                  f'ord mod P {oP} vs tropical min {vmin}  -> {"PINNED" if ok else "not pinned"}', flush=True)
        print(f'  ({time.time()-t0:.0f}s)', flush=True)
    print(f'm={m}: (a),(b) certified exactly at both points: {allok}')
