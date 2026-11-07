"""
Step 40: pointwise comparison for the angular bound.
lam_ours = min Rayleigh quotient of  -(p S')' + U S,  p = (1+al x^2)(1-x^2),  U = -V (V as in step37).
Split U = U_std + Sing_al + E with
  U_std  = (m^2 + s^2 + 2 s m x)/(1-x^2) + s - s^2          (a = 0 spin-weighted operator, lowest = 2 for s = +-2)
  Sing_al = al [m^2 (2+al-x^2)/(1+al x^2) + s^2 + 2 s m x]/(1-x^2)  >= 0
  E = -s al + 2 al x^2 - s^2 al (1-x^2)/(1+al x^2) + (1+al)^2/(1+al x^2) [xi^2 (1-x^2) + 2 xi (m + s x)]
      - 2 s al m (1+al) x/(1+al x^2)
Since p >= 1-x^2 and Sing_al >= 0:  lam_ours >= 2 + inf_x E.   Survey inf E on a superset of the window:
xi = a omega, omega/m in [0, Omega_H], Omega_H = a/(r1^2 + a^2).
"""
import numpy as np
from step3_dual import kds_roots
from step9_freeq import valid_kds

def E(x, al, xi, m, s):
    return (-s*al + 2*al*x**2 - s**2*al*(1 - x**2)/(1 + al*x**2)
            + (1 + al)**2/(1 + al*x**2)*(xi**2*(1 - x**2) + 2*xi*(m + s*x))
            - 2*s*al*m*(1 + al)*x/(1 + al*x**2))

def check_split(al, xi, m, s):   # verify U = U_std + Sing + E against step37's V at random x
    x = np.linspace(-0.97, 0.97, 7)
    V = (-s*(1 - al) + ((1 + al)**2/al)*xi**2 - 2*al*x**2
         + (1 + al)/(1 + al*x**2)*(2*s*(al*m - (1 + al)*xi)*x - ((1 + al)**2/al)*xi**2 - 2*m*(1 + al)*xi + s**2)
         - (1 + al)**2*m**2/((1 + al*x**2)*(1 - x**2)) - (1 + al)*(s**2 + 2*s*m*x)/(1 - x**2))
    Ustd = (m*m + s*s + 2*s*m*x)/(1 - x**2) + s - s*s
    Sing = al*(m*m*(2 + al - x**2)/(1 + al*x**2) + s*s + 2*s*m*x)/(1 - x**2)
    return np.max(np.abs(-V - (Ustd + Sing + E(x, al, xi, m, s))))

if __name__ == '__main__':
    print('split identity error:', max(check_split(0.03, 0.2, m, s) for m in (-3, -1, 1, 2) for s in (-2, 2)))
    xs = np.linspace(-1, 1, 2001)
    worst = {}
    for Lam in np.geomspace(1e-5, 0.11, 40):
        for a in np.linspace(0.01, 1.2, 120):
            if not valid_kds(Lam, a): continue
            r3, r0, r1, r2 = [float(q) for q in kds_roots(Lam, a, 1.0)]
            al = Lam*a*a/3; OmH = a/(r1**2 + a**2)
            for m in (1, 2, 3, 4, 8):
                for u in np.linspace(0, OmH, 25):
                    for s in (-2, 2):
                        for sg in (1, -1):        # m and omega of the same sign
                            v = E(xs, al, a*sg*m*u, sg*m, s).min()
                            if v < worst.get(m, (9,))[0]: worst[m] = (v, Lam, a, u/OmH, s, sg)
    for m in sorted(worst): print('|m| =', m, ' inf E =', worst[m])
