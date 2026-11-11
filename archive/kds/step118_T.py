# Mobius-transformed coefficients T_j of N(S), S = x0*v/(1+v); exactly one sign change => exactly one root (Descartes).
import sympy as sp, pickle
t,z0,A,B,v=sp.symbols('t z0 A B v',positive=True)
x0=z0-1; mu0=z0*x0*(t-z0); mup0=-3*z0**2+2*(1+t)*z0-t
c=[6*mu0**2, -6*mu0*mup0,
   4*(z0-1)*(z0-t)*(A+B*z0)+6*t**2*z0**2-6*t**2*z0+t**2-8*t*z0**3+12*t*z0**2-6*t*z0+3*z0**4-8*z0**3+6*z0**2,
   2*(A*(t+1-2*z0)+B*mup0+3*z0*(-3*t*z0+2*t+4*z0**2-3*z0)), 0, 2*(B+4*t-15*z0+4), 15]
Pv=sp.Poly(sp.expand(sum(c[k]*x0**k*v**k*(1+v)**(6-k) for k in range(7))),v)
Tj=[sp.expand(Pv.coeff_monomial(v**j)) for j in range(7)]
for j,T in enumerate(Tj):
    PA=sp.Poly(T,A,B); print('T%d: coefA'%j, sp.factor(PA.coeff_monomial(A)),' | coefB',sp.factor(PA.coeff_monomial(B)))
pickle.dump(Tj,open('/home/claude/kds/T118.pkl','wb'))
