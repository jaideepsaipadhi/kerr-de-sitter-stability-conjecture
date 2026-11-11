"""t -> inf OUTER limit (z = t zeta, A = t^2 alpha, B = t b) of the repaired R3 multiplier (delta = 0.7, beta = 2).
   G = psi (alpha/zeta^2 + 3)/Lam,  Lam = alpha/zeta + b - 3 zeta,  psi = 1.7 (1 - exp(-2K)),
   K = (alpha/zeta0 + b - 3 zeta0)(zeta0 - zeta)^2/(zeta0^2 (1 - zeta0)).
   C_inf = (z0-z)G m - 2m - (z0-z) m',  m = zeta^2 (1-zeta), m' = 2 zeta - 3 zeta^2
   D_inf = (1/2)(m' G + m (G' + G^2)) - (z0-z)(alpha/zeta^2 + 3)(psi - 1)
   normalised by the sums of absolute values of the terms, as in check()."""
import numpy as np
def lim(z0, al, b, n=4000):
    s = np.concatenate([np.geomspace(1e-9*z0, 1e-2*z0, n), np.linspace(1e-2*z0, z0*(1 - 1e-6), n)])
    ze = z0 - s
    Lam = al/ze + b - 3*ze; Lp = -(al/ze**2 + 3)               # d/dzeta
    K0 = (al/z0 + b - 3*z0)/(z0**2*(1 - z0)); E = np.exp(-2*K0*s**2)
    psi = 1.7*(1 - E); psip = -1.7*E*4*K0*s                    # d/dzeta (s = z0 - zeta)
    g = -Lp/Lam; gp = (-2*al/ze**3)/Lam + Lp**2/Lam**2               # d/dzeta of (al/ze^2+3)/Lam, Lam' = Lp
    G = psi*g; Gp = psip*g + psi*gp
    m = ze**2*(1 - ze); mp = 2*ze - 3*ze**2
    Ct = [s*G*m, -2*m, -s*mp]; Dt = [0.5*mp*G, 0.5*m*(Gp + G**2), s*(-Lp)*(1 - psi)]
    C = sum(Ct)/sum(abs(x) for x in Ct); D = sum(Dt)/sum(abs(x) for x in Dt)
    return C.max(), D.max(), ze[C.argmax()], ze[D.argmax()]
worst = [-9, -9, None, None]
for z0 in (0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.95, 0.99, 0.999):
    for b in (5.5, 6, 8, 12, 30, 100, 1e3, 1e5):
        if b - 3*z0 <= 0: continue
        for al in (0, 1e-3, 1e-1, 1, 10, 1e3, 1e5):
            c, d, zc, zd = lim(z0, al, b)
            if c > worst[0]: worst[0], worst[2] = c, (z0, b, al, zc)
            if d > worst[1]: worst[1], worst[3] = d, (z0, b, al, zd)
print('outer limit: worst C', round(worst[0], 4), 'at (zeta0, b, alpha, zeta) =', worst[2])
print('outer limit: worst D', round(worst[1], 4), 'at', worst[3])
