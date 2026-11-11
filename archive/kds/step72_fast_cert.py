"""Fast mod-P version of the K = 7 descent certificate (item 85), for probing geometries.  Mode builds in nmod arithmetic."""
import sys, time, numpy as np, flint
from fractions import Fraction as Fr
import step52p_modp_modes as SP
import step67_4d_from_modes as F
from step66_biv_J3 import EV, VI
P = SP.PMOD
def mode_coeffs_p(Lam, a, r0, x0, K):
    bs = [SP.run(Lam, a, Fr(1), r0, x0, (Fr(s), Fr(0)), u, K) for s, u in EV]
    VIp = [[c.numerator % P*pow(c.denominator % P, P - 2, P) % P for c in row] for row in VI]
    IU = F.IU
    def co(kr, ki):
        A = bs[0][kr]; nr, nc = A.nrows(), A.ncols(); out = [dict() for _ in range(6)]
        Ar = [np.array([[int(b[kr][i, j]) for j in range(nc)] for i in range(nr)], dtype=object) for b in bs]
        Ai = [np.array([[int(b[ki][i, j]) for j in range(nc)] for i in range(nr)], dtype=object) for b in bs]
        for k in range(6):
            Mk = sum(VIp[k][q]*((Ar[q] + IU*Ai[q]) % P) for q in range(6)) % P
            nz = np.argwhere(Mk != 0)
            for i, j in nz: out[k][(int(i), int(j))] = int(Mk[i, j])
        return out, nr, nc
    C, nU, nE = co('CR', 'CI'); G, nX, _ = co('GR', 'GI'); b0 = bs[0]
    return dict(C=C, G=G, nU=nU, nE=nE, nX=nX, nk=b0['nk'], comps=b0['comps'], monsK=b0['monsK'], mons=b0['mons'])
def cert(Lam, a, r0, x0, K=7, J=4):
    t0 = time.time(); D = mode_coeffs_p(Lam, a, r0, x0, K); t1 = time.time()
    Cm, Gm, tu, tx = F.assemble(D, K)
    low = [n for n, d in enumerate(tu) if d <= J]; high = [n for n, d in enumerate(tu) if d > J]
    xh = [n for n, d in enumerate(tx) if d >= J + 2]
    d = len(low) - F.frank(Cm) + len(high) - F.frank(Gm[np.ix_(xh, high)]) - F.frank(Gm[:, low])
    return d, t1 - t0, time.time() - t1
if __name__ == '__main__':
    d, tb, ta = cert(Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5))
    print(f'K=7 J=4 structural bound d_4 <= {d}   (build {tb:.0f}s, assembly+ranks {ta:.0f}s)')
