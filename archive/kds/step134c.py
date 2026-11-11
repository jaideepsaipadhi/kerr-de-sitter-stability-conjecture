import sympy as sp, numpy as np
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
N=sp.expand(nu.subs(z,z0-S)); f=sp.lambdify((t,z0,A,B,S),N)
rng=np.random.default_rng(1); cnt={}; ex=[]
n=0
while n<20000:
    tt=1+10**rng.uniform(-3,4); zz=1+(tt-1)*rng.uniform(0.001,0.999); x0=zz-1
    AA=10**rng.uniform(-3,6); BB=3*zz-AA/zz+10**rng.uniform(-4,np.log10(tt)+3)
    L1=AA+BB-3
    if not (2*x0*L1/(tt-1)>1 and 2*x0**2*(AA/zz+3)>=tt-1): continue
    n+=1
    Sg=x0*np.linspace(0,1,4001)[1:-1]; vals=f(tt,zz,AA,BB,Sg)
    ch=int(np.sum(np.sign(vals[1:])!=np.sign(vals[:-1]))); cnt[ch]=cnt.get(ch,0)+1
    if ch!=1 and len(ex)<5: ex.append((tt,zz,AA,BB,BB-(5.5*tt-1.25)))
print(cnt); print(ex)
