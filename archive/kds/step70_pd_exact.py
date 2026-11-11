"""EXACT independence of the Lambda-PD tangents modulo gauge at the 3-jet level, BL chart (t, r, x, phi), point (3, 3/5),
Lam = 1/20, a = 3/5, M = 1  — the descent input of item 83.
Carter/PD chart (tau, r, p, sig):  tau = (t - a phi)/Xi,  sig = -phi/(a Xi),  p = a x   (constant linear map; derived by matching).
PD metric (Griffiths–Podolsky form):  Omega = 1 - alpha p r,  S = r^2 + p^2,
 g = Omega^-2 [ -Q/S (dtau - p^2 dsig)^2 + S/Q dr^2 + S/P dp^2 + P/S (dtau + r^2 dsig)^2 ],
 P = k + 2 n p - eps p^2 + 2 alpha m p^3 - (alpha^2 k + Lam/3) p^4,  Q = k - 2 m r + eps r^2 - 2 alpha n r^3 - (alpha^2 k + Lam/3) r^4.
Background m = M, n = 0, eps = 1 - Lam a^2/3, k = a^2, alpha = 0.  Tangents h_i = pullback of dg/d(param_i) (exact rational).
Checks: (0) pulled-back background equals the BL metric of step52 exactly; then rank_Q[G_3; PD_3] - rank_Q G_3."""
import itertools, time, sympy as sp, flint
from fractions import Fraction as Fr
import step52_exact_modes as S
from step67_4d_from_modes import EXP
t0 = time.time()
r, x = sp.symbols('r x'); m_, n_, e_, k_, al_ = sp.symbols('m n eps k alpha')
Lam, a, M = sp.Rational(1, 20), sp.Rational(3, 5), sp.Integer(1); Xi = 1 + Lam*a**2/3
r0, x0 = sp.Integer(3), sp.Rational(3, 5)
p = a*x
P = k_ + 2*n_*p - e_*p**2 + 2*al_*m_*p**3 - (al_**2*k_ + Lam/3)*p**4
Q = k_ - 2*m_*r + e_*r**2 - 2*al_*n_*r**3 - (al_**2*k_ + Lam/3)*r**4
Sm = r**2 + p**2; Om = 1 - al_*p*r
# Carter-chart metric components in coordinates (tau, r, p, sig)
gC = sp.zeros(4)
gC[0, 0] = (-Q + P)/Sm; gC[0, 3] = gC[3, 0] = (Q*p**2 + P*r**2)/Sm; gC[3, 3] = (-Q*p**4 + P*r**4)/Sm
gC[1, 1] = Sm/Q; gC[2, 2] = Sm/P
gC = gC/Om**2
# Jacobian d(tau, r, p, sig)/d(t, r, x, phi)
Jm = sp.Matrix([[1/Xi, 0, 0, -a/Xi], [0, 1, 0, 0], [0, 0, a, 0], [0, 0, 0, -1/(a*Xi)]])
bg = {m_: M, n_: 0, e_: 1 - Lam*a**2/3, k_: a**2, al_: 0}
gBL = (Jm.T*gC*Jm).subs(bg).applyfunc(sp.cancel)
# compare with BL metric (as in step52): rho2 = r^2 + a^2 x^2
Dr = (r**2 + a**2)*(1 - Lam*r**2/3) - 2*M*r; Dx = 1 + Lam*a**2*x**2/3; rho2 = r**2 + a**2*x**2; s2 = 1 - x**2
ref = sp.zeros(4)
ref[0, 0] = -(Dr - Dx*a**2*s2)/(rho2*Xi**2); ref[0, 3] = ref[3, 0] = a*s2*(Dr - Dx*(r**2 + a**2))/(rho2*Xi**2)
ref[3, 3] = s2*(Dx*(r**2 + a**2)**2 - Dr*a**2*s2)/(rho2*Xi**2); ref[1, 1] = rho2/Dr; ref[2, 2] = rho2/(Dx*s2)
diff = (gBL - ref).applyfunc(sp.cancel)
print('(0) pulled-back PD background == BL KdS metric exactly:', diff == sp.zeros(4), f'[{time.time()-t0:.0f}s]', flush=True)
# exact 3-jets of the 5 tangents at (r0, x0), in the mode builder's unknown ordering
K = 4; ref_run = S.run(Fr(1, 20), Fr(3, 5), Fr(1), Fr(3), Fr(3, 5), (Fr(1), Fr(0)), 1, K)
mons, monsK, nk, comps = ref_run['mons'], ref_run['monsK'], ref_run['nk'], ref_run['comps']
J3 = [(c, j) for c in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= 3]
u, v = sp.symbols('u v')
def jets(expr):
    # exact Taylor coefficients to total degree 3 via derivatives at the point
    f = sp.cancel(expr); out = {}
    for d in range(4):
        for i in range(d + 1):
            jdx = d - i
            dv = sp.diff(f, r, i, x, jdx) if (i or jdx) else f
            out[(i, jdx)] = sp.Rational(dv.subs({r: r0, x: x0}))/(sp.factorial(i)*sp.factorial(jdx))
    return out
PDv = []
for prm in (m_, n_, e_, k_, al_):
    h = (Jm.T*sp.diff(gC, prm)*Jm).subs(bg)
    vec = []; cache = {}
    for (c, j) in J3:
        pq = comps[c]; mon = mons[monsK[j]]
        if pq not in cache: cache[pq] = jets(h[pq[0], pq[1]])
        vec.append(sp.Rational(cache[pq].get(tuple(mon), 0)))
    PDv.append(vec)
print('PD tangent 3-jets computed', f'[{time.time()-t0:.0f}s]', flush=True)
# exact gauge 3-jets in 4D from exact mode gauge coefficients (degree <= 1 in (sigma, m)): sigma <-> i d_t, m <-> -i d_phi.
# A 4D gauge generator t^j phi^k X_mu with j + k >= 1 contributes to t-independent 3-jets only through d_t / d_phi of order 1.
from step66_biv_J3 import EV, VI
bs = [S.run(Fr(1, 20), Fr(3, 5), Fr(1), Fr(3), Fr(3, 5), (Fr(s), Fr(0)), uu, K) for s, uu in EV]
nX, ncol = bs[0]['GR'].nrows(), bs[0]['GR'].ncols()
Gco = [dict() for _ in range(6)]
for i in range(nX):
    for j in range(ncol):
        vr = [Fr(int(b['GR'][i, j].p), int(b['GR'][i, j].q)) for b in bs]; vi = [Fr(int(b['GI'][i, j].p), int(b['GI'][i, j].q)) for b in bs]
        if not any(vr) and not any(vi): continue
        for kk in range(6):
            cr = sum(VI[kk][q]*vr[q] for q in range(6)); ci = sum(VI[kk][q]*vi[q] for q in range(6))
            if cr or ci: Gco[kk][(i, j)] = (cr, ci)
from math import factorial
deg = lambda i: sum(mons[i])
udeg = [deg(monsK[j]) for c in range(len(comps)) for j in range(nk)]
Xmons = [i for i in range(len(mons)) if deg(i) <= K + 1]; xdeg = [deg(Xmons[j]) for f in range(nX//len(Xmons)) for j in range(len(Xmons))]
tp = [(j, k) for j in range(5) for k in range(5 - j)]
U3 = [(j, k, i) for (j, k) in tp for i in range(ncol) if j + k + udeg[i] <= 3]; Ui = {w: n for n, w in enumerate(U3)}
X4 = [(j, k, i) for (j, k) in tp for i in range(nX) if j + k + xdeg[i] <= 4]
IP = [(1, 0), (0, 1), (-1, 0), (0, -1)]
rows = []
for (j, k, xi) in X4:
    row = {}
    for kk, (ea, eb) in enumerate(EXP):
        if ea > j or eb > k: continue
        ph = IP[ea % 4]; ph2 = IP[(3*eb) % 4]; ph = (ph[0]*ph2[0] - ph[1]*ph2[1], ph[0]*ph2[1] + ph[1]*ph2[0])
        f = (factorial(j)//factorial(j - ea))*(factorial(k)//factorial(k - eb))
        for (ii, jj), (cr, ci) in Gco[kk].items():
            if ii != xi: continue
            w = (j - ea, k - eb, jj)
            if w not in Ui: continue
            o = row.get(Ui[w], (Fr(0), Fr(0)))
            row[Ui[w]] = (o[0] + f*(ph[0]*cr - ph[1]*ci), o[1] + f*(ph[0]*ci + ph[1]*cr))
    if row: rows.append([row.get(n, (Fr(0), Fr(0))) for n in range(len(U3))])
print('4D gauge 3-jet generators', len(rows), 'x', len(U3), flush=True)
def realified(R):
    nr, nc = len(R), len(R[0]); A = flint.fmpq_mat(2*nr, 2*nc)
    for i, row in enumerate(R):
        for j, (xr, yi) in enumerate(row):
            X = flint.fmpq(xr.numerator, xr.denominator); Y = flint.fmpq(yi.numerator, yi.denominator)
            A[i, j] = X; A[nr + i, nc + j] = X; A[i, nc + j] = Y; A[nr + i, j] = -Y
    return A
rG = realified(rows).rank()//2
PDc = []
for vec in PDv:
    full = [(Fr(0), Fr(0))]*len(U3)
    for (c, j), cv in zip(J3, vec): full[Ui[(0, 0, c*nk + j)]] = (Fr(int(cv.p), int(cv.q)), Fr(0))
    PDc.append(full)
rGP = realified(rows + PDc).rank()//2
print(f'gauge 3-jet rank (exact) = {rG}; with the 5 PD tangents = {rGP}  =>  PD rank mod gauge = {rGP - rG} (need 4)  [{time.time()-t0:.0f}s]')
