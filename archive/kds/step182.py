import sympy as sp, pickle
z,t,z0,A,B,Q2,K=sp.symbols('z t z0 A B Q2 Khat',real=True)   # Q2=q1^2, K=kappa*q1
mu=z*(z-1)*(t-z); L=A/z+B-3*z; Lp=sp.diff(L,z)
corr=(z-z0)**2*(2*K*Q2+K**2*(t-z))/(z*(z-1))
E=sp.together((z-z0)*(Q2+K*(t-z))*(Q2*Lp+sp.diff(corr,z))-K*Q2*sp.diff(mu,z))
n,d=sp.fraction(E); n=sp.expand(n)
print('den',sp.factor(d)); P=sp.Poly(n,A,B,Q2,K); print('monomials (A,B,Q2,K):',sorted(P.monoms()))
C=sp.expand(sp.diff(mu,z)*(z-z0)*(Q2+K*(t-z))-2*mu*(Q2+K*(t-z)-K*(z-z0)))  # C~ * q1
print('C*q1 monomials',sorted(sp.Poly(C,Q2,K).monoms()))
pickle.dump((n,d,C),open('P182.pkl','wb'))
