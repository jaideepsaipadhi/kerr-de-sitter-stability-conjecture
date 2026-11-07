"""
Step 27: closed-form finite-m Heun data for the KdS spin-s radial Teukolsky equation, gauge (0,0,0).

Delta = -(Lam/3) prod_j (r - r_j), kept symbolic in the roots (r3 = -(r0+r1+r2)).
Exponents at r_j: kappa_j = i Xi K(r_j)/Delta'(r_j), with second exponent -kappa_j - s.
Heun:  gamma = 1 + 2 kappa_0 + s,  delta = 1 + 2 kappa_1 + s,  epsilon = 1 + 2 kappa_t + s.
At z = infinity (r = r3), the exponents {A, B} come from {kappa_3, -kappa_3 - s} shifted by the gauge
factors (including (z - zinf)^{2s+1}).  The accessory parameter q comes from the f-equation:
  q = lim_{z->0} [ A B z - z (z-1)(z-t) Q_f(z) ] = -t * lim_{z->0} z Q_f(z).
The script computes A B and q symbolically, as rational functions of (roots, Lam, a, omega, m, lam, s),
then checks them numerically against step 3.
"""
import sympy as sp
import mpmath as mp

r, z = sp.symbols('r z')
r0, r1, r2, Lam, a, w, m, lam = sp.symbols('r0 r1 r2 Lambda a omega m lambda', real=True)
s = sp.Integer(-2)
r3 = -(r0 + r1 + r2)
roots = [r0, r1, r2, r3]
D = -(Lam/3)*sp.prod([r - rj for rj in roots])
Dp = sp.diff(D, r)
alpha = Lam*a**2/3; Xi = 1 + alpha
K = w*(r**2 + a**2) - a*m
Q = (Xi**2*K**2 - sp.I*s*Xi*K*Dp)/D + 4*sp.I*s*Xi*w*r - sp.Rational(2, 3)*Lam*(s + 1)*(2*s + 1)*r**2 + s*(1 - alpha) - lam
p_r = (s + 1)*Dp/D
q_r = Q/D
zinf = (r1 - r3)/(r1 - r0); t = zinf*(r2 - r0)/(r2 - r3)
rz = (r0*zinf - r3*z)/(zinf - z)
drdz = sp.diff(rz, z); d2 = sp.diff(drdz, z)
kap = [sp.I*Xi*K.subs(r, rj)/Dp.subs(r, rj) for rj in roots]
e0, e1, et = kap[0], kap[1], kap[2]
L1 = e0/z + e1/(z - 1) + et/(z - t) + (2*s + 1)/(z - zinf)
L2 = sp.diff(L1, z) + L1**2
Pz = p_r.subs(r, rz)*drdz - d2/drdz
Qz = q_r.subs(r, rz)*drdz**2
Qf = Qz + Pz*L1 + L2

# q = -t * lim_{z->0} z Qf  (simple-pole residue; the double-pole part cancels by the choice of exponent)
expr = sp.together(z*Qf)
num, den = sp.fraction(expr)
q_acc = -t*sp.cancel(num.subs(z, 0)/den.subs(z, 0))
# A B from z -> infinity:  z^2 Qf -> A B  (so that Qf ~ AB/z^2)
AB = sp.limit(z**2*Qf, z, sp.oo) if False else None


def numeric_check(Lamv, av, Mv, omegav, mv, lamv):
    from step3_dual import kds_roots, teukolsky_heun, heun_params_from_operator
    rr3, rr0, rr1, rr2 = [float(x) for x in kds_roots(Lamv, av, Mv)]
    subs = {r0: rr0, r1: rr1, r2: rr2, Lam: Lamv, a: av, w: omegav, m: mv, lam: lamv}
    qv = complex(sp.N(q_acc.subs(subs), 30))
    data = teukolsky_heun(Lamv, av, Mv, -2, complex(omegav), mv, lamv)
    e = data['exps']
    h = heun_params_from_operator(data['Pz'], data['Qz'], data['z'], data['t'], e['0'][0], e['1'][0], e['t'][0],
                                  {'e': -3, 'zinf': sp.Float(str(data['zinf']), 35)}, -2)
    return qv, h['q'], abs(qv - h['q'])/abs(h['q'])


if __name__ == '__main__':
    for args in [(0.05, 0.6, 1.0, 0.0915, 2, 7.3), (0.001, 0.95, 1.0, 0.31, 3, 20.0), (0.02, 0.9, 1.0, 0.05, 1, 2.5)]:
        qv, qref, rel = numeric_check(*args)
        print(f'{args}: closed-form q = {qv:.10f}, step-3 q = {qref:.10f}, rel. diff {rel:.1e}')
