import sympy as sp, numpy as np
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
N=sp.expand(nu.subs(z,z0-S)); f=sp.lambdify((t,z0,A,B,S),N)
bad=0;tot=0
rng=np.random.default_rng(0)
for _ in range(20000):
    tt=1+10**rng.uniform(-3,3); zz=1+(tt-1)*rng.uniform(0.001,0.999); x0=zz-1
    AA=10**rng.uniform(-3,3)*rng.integers(0,2); BB=3*zz-AA/zz+10**rng.uniform(-4,3)
    Sg=np.linspace(0,x0,4001)[1:-1]; vals=f(tt,zz,AA,BB,Sg)
    ch=np.sum(np.sign(vals[1:])!=np.sign(vals[:-1])); tot+=1
    if ch!=1: bad+=1; 
    if ch>=2 and bad<400: print("MULTI",tt,zz,AA,BB,ch)
print(bad,tot)
