"""Step 45: are the Plebanski-Demianski tangents the whole 4-dim space of step 44?
PD (vacuum, Lam) in (tau, r, p, sigma):  Omega = 1 - alpha p r,
 ds^2 = Omega^-2 [ -Q/(r^2+p^2) (dtau - p^2 dsig)^2 + (r^2+p^2)/Q dr^2 + (r^2+p^2)/P dp^2 + P/(r^2+p^2) (dtau + r^2 dsig)^2 ]
 P = k + 2 n p - eps p^2 + 2 alpha m p^3 - (alpha^2 k + Lam/3) p^4
 Q = k - 2 m r + eps r^2 - 2 alpha n r^3 - (alpha^2 k + Lam/3) r^4
Background n = alpha = 0 (Carter form of KdS).  h_i = d g / d(m, n, eps, k, alpha) by central differences.
Report rank of span{pi_J h_i} modulo gauge (expect 4 = d_J) and the constraint residuals of h_i."""
import sys, numpy as np
from step44_jet_completeness import run, counts

Lam, mP, epsP, kP = 0.05, 1.0, 1.0 - 0.05*0.36/3*0 + 0.0, 0.36   # k ~ a^2; eps chosen so that P, Q > 0 at the point
r0, p0 = 3.0, 0.3
def pd(E, mulb, recip, prm):
    m, n, eps, k, al = prm; M = E['M']
    r = E['var'](1, r0); p = E['var'](2, p0); one = E['const'](1)
    rr, pp = mulb(r, r), mulb(p, p)
    P = k*one + 2*n*p - eps*pp + 2*al*m*mulb(pp, p) - (al*al*k + Lam/3)*mulb(pp, pp)
    Q = k*one - 2*m*r + eps*rr - 2*al*n*mulb(rr, r) - (al*al*k + Lam/3)*mulb(rr, rr)
    S = rr + pp; iS = recip(S); Om = one - al*mulb(p, r); iO2 = recip(mulb(Om, Om))
    A = -mulb(Q, iS); B = mulb(P, iS)
    g = [[np.zeros(M) for _ in range(4)] for _ in range(4)]
    g[0][0] = A + B
    g[0][3] = g[3][0] = -mulb(A, pp) + mulb(B, rr)
    g[3][3] = mulb(A, mulb(pp, pp)) + mulb(B, mulb(rr, rr))
    g[1][1] = mulb(S, recip(Q)); g[2][2] = mulb(S, recip(P))
    g = [[mulb(iO2, g[i][j]) for j in range(4)] for i in range(4)]
    z = np.zeros(M)
    lv = [rr, Q, z, -one]; nv = [rr, -Q, z, -one]
    mv = [pp.astype(complex), z.astype(complex), 1j*P, one.astype(complex)]
    return g, lv, nv, mv, P[0], Q[0]

base = (mP, 0.0, epsP, kP, 0.0)
def metric_fn(E, mulb, recip):
    g, lv, nv, mv, P0, Q0 = pd(E, mulb, recip, base); print(f'   P(p0)={P0:.4f} Q(r0)={Q0:.4f} (need >0)'); return g, lv, nv, mv
def extra(E, mulb, recip):
    out = []
    for i in range(5):
        d = 1e-5
        hp = list(base); hm = list(base); hp[i] += d; hm[i] -= d
        gp = pd(E, mulb, recip, hp)[0]; gm = pd(E, mulb, recip, hm)[0]
        out.append([[(gp[a][b] - gm[a][b])/(2*d) for b in range(4)] for a in range(4)])
    return out

K = int(sys.argv[1]) if len(sys.argv) > 1 else 6
res = run(Lam, 0, 0, K, 0, 0, metric_fn=metric_fn, extra=extra)
print(f'checks: tetrad {res["chk"]:.1e}, Ric-Lam g {res["bric"]:.1e}, bg psi0 {res["bpsi0"]:.1e}, gauge {res["gauge_resid"]:.1e}')
H = res['extras']; C = res['Cmat']
print('constraint residual of the 5 PD tangents (relative):', [f'{np.abs(h @ C).max()/np.abs(h).max()/np.abs(C).max():.1e}' for h in H])
out, _ = counts(res, K)
print('d_J on this chart:', out)
monsK, mons, nk, ncomp = res['monsK'], res['E']['mons'], len(res['monsK']), len(res['comps'])
def rk(Mx, tol=1e-8):
    s_ = np.linalg.svd(Mx, compute_uv=False); return int((s_ > tol*s_[0]).sum())
for J in range(3, K - 1):
    sel = [k*nk + j for k in range(ncomp) for j, i in enumerate(monsK) if sum(mons[i]) <= J]
    Gs = res['Gh'][:, sel]; Hs = H[:, sel]
    print(f'J={J}: rank(PD tangents mod gauge) =', rk(np.vstack([Gs, Hs])) - rk(Gs))
