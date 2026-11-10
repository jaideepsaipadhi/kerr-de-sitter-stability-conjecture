import sympy as sp, pickle
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s_=z0-z; mus=z*(z-1)*(t-z); Ls=A/z+B-3*z; F=mus/s_**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mus*sp.diff(F,z),z)-2*s_*sp.diff(Ls*F,z))))
Bm=3*z0-A/(2*z0)+2*t
g=sp.expand(sp.expand(nu.subs(z,z0-S)).subs(B,Bm)*2*z0)
g1=sp.factor(g.coeff(A)); g0=sp.factor(g.subs(A,0))
zh=z0-S; muh=zh*(zh-1)*(t-zh); H=2*muh+S*sp.diff(mus,z).subs(z,zh)
Q=sp.together((H/S)*(A/z+Bm-3*z)-2*mus*(A/z**2+3))
Qn,Qd=sp.fraction(Q); Qn=sp.expand(Qn)
q1=sp.factor(Qn.coeff(A)); q0=sp.factor(Qn.subs(A,0))
print('g1',g1); print('g0',g0); print('Qden',sp.factor(Qd)); print('q1',q1); print('q0',q0)
pickle.dump((g0,g1,q0,q1,Qd),open('P172.pkl','wb'))
