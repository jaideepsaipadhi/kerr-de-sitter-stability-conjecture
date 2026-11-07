"""Decisive test for the sigma = 0 anomaly.  Take the exact mode solution space S (sigma = 0, m = 1) from step52 and lift each
basis element to a 4-D jet h = e^{i phi} H(r, cos th) in (t, r, th, phi) coordinates.  Then test it against the INDEPENDENT 4-D
constraint matrix of step44 (float).  Lifted elements violating the 4-D equations => the mode code is under-constrained (bug)."""
import numpy as np, sympy as sp
from fractions import Fraction as Fr
import flint
import step52_exact_modes as S
from step44_jet_completeness import run as run4, build
K = 4
x0 = Fr(3, 5); th0 = float(np.arccos(0.6))
import sys
SIG = (Fr(sys.argv[1]), Fr(sys.argv[2])) if len(sys.argv) > 2 else (Fr(0), Fr(0))
res = S.run(Fr(1, 20), Fr(3, 5), Fr(1), Fr(3), x0, SIG, 1, K)
CR, CI = res['CR'], res['CI']; ncol = res['ncol']
R = S.realify(CR, CI)                    # rows: [Re c ; Im c], cols: [Re eq ; Im eq]
# exact null space of the realified constraints: left kernel  v R = 0  <=>  R^T v^T = 0
RT = R.transpose()
Rr, rk = RT.rref()
ncols = RT.ncols(); piv = []
row = 0
for j in range(ncols):
    if row < rk and Rr[row, j] != 0: piv.append(j); row += 1
free = [j for j in range(ncols) if j not in piv]
basis = []
for f in free:
    v = [0.0]*ncols; v[f] = 1.0
    for i, pj in enumerate(piv): v[pj] = -float(Rr[i, f])
    basis.append(v)
nul = len(basis)
print('dim S (real) =', nul)
# 4-D constraint matrix at the same point
r4 = run4(0.05, 0.6, 1.0, K, 3.0, th0)
C4 = r4['Cmat']; E4 = r4['E']; mons4 = E4['mons']; idx4 = E4['idx']; monsK4 = r4['monsK']; comps4 = r4['comps']
mulb = E4['mulb']; N4 = K + 1
s_, c_ = E4['sincos'](E4['var'](2, th0)); yx = c_ - float(x0)*E4['const'](1)
yr = E4['var'](1, 3.0) - 3.0*E4['const'](1)
ang = E4['var'](3, 0.0) - float(SIG[0])*(E4['var'](0, 0.0))
sphi, cphi = E4['sincos'](ang)
# exp(si * yt)
ex = E4['const'](1); term = E4['const'](1); yt = E4['var'](0, 0.0)
for k in range(1, K + 2):
    term = mulb(term, yt)*float(SIG[1])/k; ex = ex + term
sphi, cphi = mulb(sphi, ex), mulb(cphi, ex)
monsM, nkM, compsM = res['mons'], res['nk'], res['comps']; monsKM = res['monsK']
def pw(f, k):
    out = E4['const'](1)
    for _ in range(k): out = mulb(out, f)
    return out
def lift(vec):
    re, im = vec[:ncol], vec[ncol:]
    Hc = {}
    for k0, pq in enumerate(compsM):
        Fr_ = np.zeros(E4['M']); Fi = np.zeros(E4['M'])
        for j, mi in enumerate(monsKM):
            i_, j_ = monsM[mi]
            term = mulb(pw(yr, i_), pw(yx, j_))
            Fr_ += re[k0*nkM + j]*term; Fi += im[k0*nkM + j]*term
        Hc[pq] = (Fr_, Fi)
    # h4 (real part of e^{i phi} H): Re = Fr cos - Fi sin ; in theta basis: d x = -s d th
    out = []
    for part in (0, 1):
        h4 = {}
        for (p, q), (Fr_, Fi) in Hc.items():
            val = (mulb(Fr_, cphi) - mulb(Fi, sphi)) if part == 0 else (mulb(Fr_, sphi) + mulb(Fi, cphi))
            fac = E4['const'](1)
            for idxc in (p, q):
                if idxc == 2: fac = mulb(fac, -s_)
            h4[(p, q)] = mulb(val, fac)
        v = np.zeros(C4.shape[0]); nk4 = len(monsK4)
        for k0, pq in enumerate(comps4):
            if pq in h4:
                v[k0*nk4:(k0 + 1)*nk4] = h4[pq][monsK4]
        out.append(v)
    return out
worst = 0
for b in basis:
    for v in lift(b):
        if np.abs(v).max() == 0: continue
        r_ = np.abs(v @ C4).max()/(np.abs(v).max()*np.abs(C4).max())
        worst = max(worst, r_)
print('max relative 4-D constraint residual over lifted mode solutions:', worst)
# control: lift a pure mode-GAUGE element, which must satisfy the 4-D constraints
g0 = [float(res['GR'][5, j]) for j in range(ncol)] + [float(res['GI'][5, j]) for j in range(ncol)]
print('control (lifted gauge row) residual:', max(np.abs(v @ C4).max()/(np.abs(v).max()*np.abs(C4).max() + 1e-300) for v in lift(g0)))
keep4 = [k for k, mo in enumerate(mons4) if sum(mo) <= K - 2]; nb = len(keep4)
names = [f'Ric{pq}' for pq in comps4] + ['psi0.re', 'psi0.im', 'psi4.re', 'psi4.im']
blk = np.zeros(len(names))
for b in basis:
    for v in lift(b):
        if np.abs(v).max() == 0: continue
        rr_ = v @ C4
        for i in range(len(names)):
            blk[i] = max(blk[i], np.abs(rr_[i*nb:(i + 1)*nb]).max()/(np.abs(v).max()*np.abs(C4).max()))
print('violation by 4-D equation block:'); [print('  ', n, f'{x:.2e}') for n, x in zip(names, blk)]
# localise Ricci violations by 4-D monomial type (which of t, r, th, phi appear)
from collections import defaultdict
byt = defaultdict(float)
for b in basis[:60]:
    for v in lift(b):
        if np.abs(v).max() == 0: continue
        rr_ = v @ C4
        for i in range(10):
            for jj, kk in enumerate(keep4):
                mo = mons4[kk]; key = ''.join(n for n, e in zip('t r th ph'.split(), mo) if e) or 'const'
                byt[key] = max(byt[key], abs(rr_[i*nb + jj])/(np.abs(v).max()*np.abs(C4).max()))
print('Ricci violation by monomial type:', {k: f'{x:.1e}' for k, x in sorted(byt.items())})
