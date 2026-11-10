import numpy as np, sympy as sp, mpmath as mp
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s_=z0-z; mus=z*(z-1)*(t-z); Ls=A/z+B-3*z; F=mus/s_**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mus*sp.diff(F,z),z)-2*s_*sp.diff(Ls*F,z))))
NP=sp.Poly(sp.expand(nu.subs(z,z0-S)),S); cf=[sp.lambdify((t,z0,A,B),c) for c in NP.all_coeffs()]
def pts(n,seed):
    rng=np.random.default_rng(seed); out=[]
    while len(out)<n:
        tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-3,1-1e-3); x0=zz-1
        a0=max(0.0,zz*((tt-1)/(2*x0**2)-3)); AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
        bb=AA/(2*zz)+2*zz+10**rng.uniform(-4,1); BB=3*zz-AA/zz+bb
        if 2*x0*(AA+BB-3)<=tt-1: continue
        c=np.array([g(tt,zz,AA,BB) for g in cf],float); r=np.roots(c)
        r=[x.real for x in r if abs(x.imag)<1e-9*abs(x) and 0<x.real<x0]
        if len(r)!=1: continue
        out.append((tt,zz,AA,BB,r[0]))
    return out
def cond(P,kf,dkf):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); w0=(2*mu_+Sh*mup)/Sh
    sg=np.linspace(1e-6,1,1500); q=zz-Sh*sg; s=Sh*sg; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    k=w0*kf(sg); kz=-w0*dkf(sg)/Sh
    return (mu*kz+k**2-2*s*L*k+2*mu*s*Lp).max()<=0
P=pts(9000,23)
for c in (1.0,):
  for p in (2.0,):
    kf=lambda x,p=p,c=c: x**p*(1-c*x*(1-x)); dkf=lambda x,p=p,c=c: p*x**(p-1)*(1-c*x*(1-x))+x**p*(-c*(1-2*x))
    print(c,p,sum(cond(Q,kf,dkf) for Q in P),'/',len(P),flush=True)
c=1.0; kf=lambda x: x**2*(1-c*x*(1-x)); dkf=lambda x: 2*x*(1-c*x*(1-x))+x**2*(-c*(1-2*x))
for Q in P:
    if not cond(Q,kf,dkf):
        tt,zz,AA,BB,Sh=Q; x0=zz-1; print('t',tt,'zf',x0/(tt-1),'A',AA,'beta',AA/zz+BB-3*zz,'m',AA/(2*zz)+zz,'Sh/x0',Sh/x0)
