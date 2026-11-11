import sympy as sp, numpy as np, pickle
z,t,z0,A,B,S,be=sp.symbols('z t z0 A B S beta')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
N=sp.expand(nu.subs(z,z0-S).subs(B,3*z0-A/z0+be)*z0)
cA=sp.factor(N.coeff(A)); cb=sp.factor(N.coeff(be)); print('A-coef',cA); print('beta-coef',cb)
H=(2*mu+s*sp.diff(mu,z)).subs(z,z0-S); zz=z0-S
AK=sp.together(zz**2*(H-3*S**3)/S**3)
N0=N.subs(be,0)
pickle.dump((N0,AK,cA,cb),open('P2_136.pkl','wb'))
f=sp.lambdify((t,z0,S,A),N0); fk=sp.lambdify((t,z0,S),AK); fa=sp.lambdify((t,z0,S),cA)
rng=np.random.default_rng(3); bad=0; badA=0
for _ in range(100000):
    tt=1+10**rng.uniform(-2,3); z_=1+(tt-1)*rng.uniform(1e-3,1-1e-3); x0=z_-1; SS=x0*rng.uniform(1e-3,1-1e-3)
    if fa(tt,z_,SS)>0: badA+=1
    a=max(fk(tt,z_,SS),0.0)
    if f(tt,z_,SS,a)>=0: bad+=1
print('Acoef>0',badA,'N>=0 at corner',bad)
bad=0;n=0;ex=[]
for _ in range(100000):
    tt=1+10**rng.uniform(-2,3); z_=1+(tt-1)*rng.uniform(1e-3,1-1e-3); x0=z_-1; SS=x0*rng.uniform(1e-3,1-1e-3)
    a0=max(0.0,z_*((tt-1)/(2*x0**2)-3)); a=max(fk(tt,z_,SS),a0); n+=1
    if f(tt,z_,SS,a)>=0: bad+=1; ex.append((tt,z_,SS/x0,a0,fk(tt,z_,SS)))
print('with R3 gate',n,bad,ex[:4])
