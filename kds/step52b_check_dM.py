"""Consistency test for the sigma = 0 anomaly: h = dg/dM (exact stationary axisymmetric solution, dpsi0 = dpsi4 = 0)
must satisfy c . C = 0 exactly in the mode code at sigma = 0, m = 0."""
from fractions import Fraction as Fr
import flint, step52_exact_modes as S
K = 5
Lam, a, Mm, r0, x0 = Fr(1, 20), Fr(3, 5), Fr(1), Fr(3), Fr(3, 5)
res = S.run(Lam, a, Mm, r0, x0, (Fr(0), Fr(0)), 0, K)
C = S.Ctx(K + 1)
r = C.var(0, r0); x = C.var(1, x0); one = C.const(1)
al = Lam*a*a/3; Xi = 1 + al
rr = C.mul(r, r); xx = C.mul(x, x); s2 = C.add(one, C.scal(-1, xx)); rho2 = C.add(rr, C.scal(a*a, xx))
rpa = C.add(rr, C.const(a*a))
Dr = C.add(C.mul(rpa, C.add(one, C.scal(-Lam/3, rr))), C.scal(-2*Mm, r))
A1p = C.scal(2/Xi**2, C.mul(r, C.recip(rho2)))          # dA1/dM = 2r/(Xi^2 rho^2)
h = {(0, 0): A1p, (0, 3): C.scal(-a, C.mul(A1p, s2)), (3, 3): C.scal(a*a, C.mul(A1p, C.mul(s2, s2))),
     (1, 1): C.mul(C.scal(2, C.mul(r, rho2)), C.recip(C.mul(Dr, Dr)))}
comps, nk, monsK = res['comps'], res['nk'], res['monsK']
v = flint.fmpq_mat(1, res['ncol'])
for k0, pq in enumerate(comps):
    if pq in h:
        for j, mi in enumerate(monsK): v[0, k0*nk + j] = flint.fmpq(h[pq][mi].numerator, h[pq][mi].denominator)
R = v*res['CR']; I = v*res['CI']
nz = sum(1 for j in range(R.ncols()) if R[0, j] != 0 or I[0, j] != 0)
print('dg/dM against the exact sigma=0 constraints: nonzero residual entries =', nz, 'of', R.ncols())
