import numpy as np, sympy as sp
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s_=z0-z; mus=z*(z-1)*(t-z); Ls=A/z+B-3*z; F=mus/s_**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mus*sp.diff(F,z),z)-2*s_*sp.diff(Ls*F,z))))
Nf=sp.lambdify((t,z0,A,B,S),sp.expand(nu.subs(z,z0-S)))
ths=np.linspace(0.3,1.95,34)
rng=np.random.default_rng(3); tot=0; ok=0; okstart={}
for i in range(30000):
    tt=1+10**rng.uniform(-2,3.5); zz=1+(tt-1)*10**rng.uniform(-3,-0.0001); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3)); AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=AA/(2*zz)+2*zz+10**rng.uniform(-4,2); BB=3*zz-AA/zz+bb
    if 2*x0*(AA+BB-3)<=tt-1: continue
    u=np.linspace(1e-6,1-1e-6,600); q=1+x0*u; s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    mup=(q-1)*(tt-q)+q*(tt-q)-q*(q-1); H=2*mu+s*mup; N=Nf(tt,zz,AA,BB,s)
    tot+=1; good=False
    for th in ths:
        if (th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp).min()<=0: continue
        fail=H-th*s*s*L<0
        if not np.any(fail & (N>=0)): good=True; break
    ok+=good
print(tot,ok)
