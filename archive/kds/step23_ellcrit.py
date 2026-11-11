"""
Step 23: n1 = -mu(z) [Qz2(z) + P1(z)^2/4] should be independent of z.  Write n1 = C0 (ell - ell_crit).
Check that it is constant in z, and compute ell_crit across the window.  If ell_crit <= 0 then n1 > 0
for every ell > 0.
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, mpmath as mp
from step3_dual import kds_roots
from step5_window import ctdc_window
mp.mp.dps = 40

def n1_fun(Lam, a, w, ell, zz):
    r3, r0, r1, r2 = [mp.mpf(x) for x in kds_roots(Lam, a, 1.0)]
    Lam, a, w, ell = [mp.mpf(x) for x in (Lam, a, w, ell)]
    zinf = (r1 - r3)/(r1 - r0); t = zinf*(r2 - r0)/(r2 - r3)
    Xi = 1 + Lam*a**2/3
    D = lambda r: (r**2 + a**2)*(1 - Lam*r**2/3) - 2*r
    Dp = lambda r: 2*r*(1 - Lam*r**2/3) - (r**2 + a**2)*2*Lam*r/3 - 2
    Kh = lambda r: w*(r**2 + a**2) - a
    rz = lambda z: (r0*zinf - r3*z)/(zinf - z)
    drdz = lambda z: zinf*(r0 - r3)/(zinf - z)**2
    k = [1j*Xi*Kh(rj)/Dp(rj) for rj in (r0, r1, r2)]
    P1 = lambda z: 2*k[0]/z + 2*k[1]/(z - 1) + 2*k[2]/(z - t)
    Qz2 = lambda z: (Xi**2*Kh(rz(z))**2/D(rz(z))**2 - ell/D(rz(z)))*drdz(z)**2
    mu = lambda z: z*(z - 1)*(t - z)
    return [complex(-mu(z)*(Qz2(z) + P1(z)**2/4)) for z in (mp.mpf(x)*1 for x in zz)], float(t)

for (Lam, a) in [(0.05, 0.6), (0.001, 0.95), (0.05, 0.9877), (0.001, 0.999), (0.0001, 0.99), (0.11, 0.01)]:
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    out = []
    for wf in (0.001, 0.25, 0.5, 0.75, 0.999):
        w = lo + wf*(hi - lo)
        _, t = n1_fun(Lam, a, w, 0.0, [1.5])
        zs = [1 + (t - 1)*f for f in (0.13, 0.5, 0.87)] + [-3.7, 0.4]
        v0, _ = n1_fun(Lam, a, w, 0.0, zs)
        v1, _ = n1_fun(Lam, a, w, 1.0, zs)
        spread = max(abs(x - v0[0]) for x in v0)/max(abs(v0[0]), abs(v1[0] - v0[0]))
        C0 = (v1[0] - v0[0]).real
        ell_crit = -v0[0].real/C0
        out.append(f'wf={wf}: ell_crit={ell_crit:+.3e} (C0={C0:.4f}, z-spread {spread:.0e}, Im {abs(v0[0].imag):.0e})')
    print(f'Lam={Lam} a={a}:'); [print('    ', o) for o in out]
