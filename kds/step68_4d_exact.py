"""RIGOROUS upper bound on the 4D jet count d_J:  d_J = dim pi_J(S) - dim pi_J(G),
   dim pi_J(S) = n_L - rank C + rank C_H   (C = full constraint matrix, C_H = its rows of degree > J)   -- EXACT ranks over Q(i)
   dim pi_J(G) >= rank_P G_L                                                                          -- mod-P lower bound
so d_J <= n_L - rank_Q C + rank_Q C_H - rank_P G_L.   Complex ranks via realification (rank_C = rank_Q[[Re, Im], [-Im, Re]] / 2).
Exact 4D assembly from exact mode coefficient matrices (sigma <-> i d_t, m <-> -i d_phi on t^j phi^k)."""
import sys, time, numpy as np, flint
from fractions import Fraction as Fr
from math import factorial
import step52_exact_modes as S
from step66_biv_J3 import EV, VI
import step67_4d_from_modes as F

EXP = F.EXP
def mode_coeffs_exact(Lam, a, r0, x0, K):
    bs = [S.run(Lam, a, Fr(1), r0, x0, (Fr(s), Fr(0)), u, K) for s, u in EV]
    A = bs[0]['CR']; nr, nc = A.nrows(), A.ncols(); C = [dict() for _ in range(6)]
    for i in range(nr):
        for j in range(nc):
            vr = [Fr(int(b['CR'][i, j].p), int(b['CR'][i, j].q)) for b in bs]
            vi = [Fr(int(b['CI'][i, j].p), int(b['CI'][i, j].q)) for b in bs]
            if not any(vr) and not any(vi): continue
            for k in range(6):
                cr = sum(VI[k][q]*vr[q] for q in range(6)); ci = sum(VI[k][q]*vi[q] for q in range(6))
                if cr or ci: C[k][(i, j)] = (cr, ci)
    b0 = bs[0]
    return dict(C=C, nU=nr, nE=nc, nk=b0['nk'], comps=b0['comps'], monsK=b0['monsK'], mons=b0['mons'])

IPOW = [(1, 0), (0, 1), (-1, 0), (0, -1)]          # i^k
def cmul(p, q): return (p[0]*q[0] - p[1]*q[1], p[0]*q[1] + p[1]*q[0])

def assemble_exact(D, K, J):
    mons = D['mons']; deg = lambda i: sum(mons[i]); monsK, nk, ncomp = D['monsK'], D['nk'], len(D['comps'])
    keep = [i for i in range(len(mons)) if deg(i) <= K - 2]; nblk = D['nE']//len(keep)
    udeg = [deg(monsK[j]) for c in range(ncomp) for j in range(nk)]
    edeg = [deg(keep[j]) for b in range(nblk) for j in range(len(keep))]
    tp = [(j, k) for j in range(K + 2) for k in range(K + 2 - j)]
    U = [(j, k, i) for (j, k) in tp for i in range(D['nU']) if j + k + udeg[i] <= K]
    Eq = [(j, k, e) for (j, k) in tp for e in range(D['nE']) if j + k + edeg[e] <= K - 2]
    Ui = {u: n for n, u in enumerate(U)}; Ei = {e: n for n, e in enumerate(Eq)}
    ent = {}
    for ku, (a, b) in enumerate(EXP):
        ph = cmul(IPOW[a % 4], IPOW[(3*b) % 4])           # i^a (-i)^b
        for (i, e), c in D['C'][ku].items():
            cc = cmul(ph, c)
            for (j, k) in tp:
                if a > j or b > k: continue
                u = (j, k, i); q = (j - a, k - b, e)
                if u in Ui and q in Ei:
                    f = (factorial(j)//factorial(j - a))*(factorial(k)//factorial(k - b))
                    key = (Ui[u], Ei[q]); old = ent.get(key, (Fr(0), Fr(0)))
                    ent[key] = (old[0] + f*cc[0], old[1] + f*cc[1])
    tot = [j + k + udeg[i] for (j, k, i) in U]
    return ent, len(U), len(Eq), tot

def realified(ent, rows, nE):
    ri = {r: n for n, r in enumerate(rows)}; nr = len(rows)
    M = flint.fmpq_mat(2*nr, 2*nE)
    for (u, e), (x, y) in ent.items():
        if u not in ri: continue
        a_ = ri[u]
        X = flint.fmpq(x.numerator, x.denominator); Y = flint.fmpq(y.numerator, y.denominator)
        if x: M[a_, e] = X; M[nr + a_, nE + e] = X
        if y: M[a_, nE + e] = Y; M[nr + a_, e] = -Y
    return M

if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 6; Js = [int(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [3]
    t0 = time.time(); Lam, a, r0, x0 = Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)
    D = mode_coeffs_exact(Lam, a, r0, x0, K); print('exact mode coefficients', round(time.time() - t0), 's', flush=True)
    Dm = F.mode_coeffs(Lam, a, r0, x0, K); Cm, Gm, tu, tx = F.assemble(Dm, K)
    ent, nU, nE, tot = assemble_exact(D, K, None); print('exact 4D assembly', nU, 'x', nE, round(time.time() - t0), 's', flush=True)
    rC = realified(ent, list(range(nU)), nE).rank()//2; print('exact rank C =', rC, round(time.time() - t0), 's', flush=True)
    for J in Js:
        high = [n for n, d in enumerate(tot) if d > J]; low = [n for n, d in enumerate(tot) if d <= J]
        rCH = realified(ent, high, nE).rank()//2
        rGL = F.frank(Gm[:, low])
        print(f'J={J}: n_L={len(low)}, rank C={rC}, rank C_H={rCH}, rank_P G_L={rGL}  =>  d_J <= {len(low) - rC + rCH - rGL}  (RIGOROUS)  [{round(time.time() - t0)} s]', flush=True)
