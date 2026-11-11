import sympy as sp, numpy as np
z,t,z0,A,B,g=sp.symbols('z t z0 A B gamma')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
N=sp.expand(nu)
H=2*mu+s*sp.diff(mu,z); T=sp.expand(H-s**2*L)
Bsub=A/(2*z0)+2*t+g-A/z0+3*z0     # ell=A/z0+B-3z0 = A/(2z0)+2t+gamma
Ng=sp.expand(N.subs(B,Bsub)); Tg=sp.expand(T.subs(B,Bsub))
NB=sp.factor(N.coeff(B)); print('N_B',NB)
cg_T=Tg.coeff(g); cg_N=Ng.coeff(g); print('T_gamma',sp.factor(cg_T),' N_gamma',sp.factor(cg_N))
kap=sp.factor(-cg_T/cg_N)   # minimal kappa making gamma-coef zero
C=sp.factor(sp.expand(Tg-kap*Ng))
cA=sp.factor(sp.expand(sp.numer(sp.together(C)).coeff(A))); c0=sp.factor(sp.numer(sp.together(C)).subs({A:0,g:0}))
den=sp.factor(sp.denom(sp.together(C)))
print('kappa',kap); print('den',den); print('cA',cA); print('c0',c0)
import pickle; pickle.dump((kap,C,cA,c0,den),open('P154.pkl','wb'))
