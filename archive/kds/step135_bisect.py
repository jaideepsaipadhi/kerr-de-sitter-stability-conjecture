import sympy as sp, numpy as np, pickle
z,t,z0,A,B,v,be=sp.symbols('z t z0 A B v beta')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
x0=z0-1; Bs=3*z0-A/z0+be
def coeffs(a,b):  # S in (a,b) -> S=(a+b v)/(1+v)
    N=sp.expand(sp.cancel(nu.subs(z,z0-(a+b*v)/(1+v))*(1+v)**6).subs(B,Bs)*z0)
    P=sp.Poly(N,v); return [sp.expand(P.coeff_monomial(v**k)) for k in range(7)]
TL=coeffs(0,x0/2); TR=coeffs(x0/2,x0)
pickle.dump((TL,TR),open('T135.pkl','wb'))
fL=sp.lambdify((t,z0,A,be),TL); fR=sp.lambdify((t,z0,A,be),TR)
rng=np.random.default_rng(9); cnt={}
def ch(c):
    c=np.array(c,float); sc=np.abs(c).max(); c=c[np.abs(c)>1e-12*sc]; sg=np.sign(c); return int(np.sum(sg[1:]!=sg[:-1]))
for _ in range(100000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0_=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0_**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1); need=(tt-1)/(2*x0_)-3*x0_-AA*x0_/zz
    bb=max(need,0)+10**rng.uniform(-3,3)*(1+abs(need))
    k=(ch(fL(tt,zz,AA,bb)),ch(fR(tt,zz,AA,bb))); cnt[k]=cnt.get(k,0)+1
print(cnt)
