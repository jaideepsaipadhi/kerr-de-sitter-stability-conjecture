import sympy as sp, numpy as np, mpmath as mp
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
NP=sp.Poly(sp.expand(nu.subs(z,z0-S)),S); cf=[sp.lambdify((t,z0,A,B),c) for c in NP.all_coeffs()]
mp.mp.dps=40
rng=np.random.default_rng(4); res=[]; neg=0
for i in range(3000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=10**rng.uniform(-4,2)*(1+AA/zz)  # beta = L(z0)
    BB=3*zz-AA/zz+bb
    if 2*x0*(AA+BB-3)<=tt-1: continue
    c=[mp.mpf(g(tt,zz,AA,BB)) for g in cf]
    r=[mp.re(x) for x in mp.polyroots(c,maxsteps=400,extraprec=400) if abs(mp.im(x))<1e-20 and 0<mp.re(x)<x0]
    if len(r)!=1: print('P1?',len(r)); continue
    Sh=r[0]; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1)
    H=2*mu_+Sh*mup
    Q=H/Sh-2*(AA*(zz/zh-1-mp.log(zz/zh))+1.5*Sh**2)
    res.append(float(Q/(H/Sh))); neg+=Q<=0
print(len(res),'neg',neg,'min normalised',min(res))
from scipy.integrate import solve_ivp
rng=np.random.default_rng(4); out=[]
for i in range(3000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=10**rng.uniform(-4,2)*(1+AA/zz); BB=3*zz-AA/zz+bb
    if 2*x0*(AA+BB-3)<=tt-1: continue
    c=[mp.mpf(g(tt,zz,AA,BB)) for g in cf]
    r=[mp.re(x) for x in mp.polyroots(c,maxsteps=400,extraprec=400) if abs(mp.im(x))<1e-20 and 0<mp.re(x)<x0]
    Sh=float(r[0]); zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1)
    H=2*mu_+Sh*mup; Q=H/Sh-2*(AA*(zz/zh-1-np.log(zz/zh))+1.5*Sh**2)
    if Q>0: continue
    # y=(f, P=mu f' - 2 s L f); P' = -2 L f  (since (mu f')'-2s(Lf)' => (mu f' -2sLf)' = -2Lf ... check: d/dz(2sLf)=2s(Lf)'-2Lf)
    Lf=lambda zq: AA/zq+BB-3*zq; muf=lambda zq: zq*(zq-1)*(tt-zq)
    f0=mu_/Sh**2; fp0=f0*H/(Sh*mu_)
    P0=mu_*fp0-2*Sh*Lf(zh)*f0
    def rhs(zq,y):
        f,P=y; sq=zz-zq; fp=(P+2*sq*Lf(zq)*f)/muf(zq); return [fp,2*Lf(zq)*f]
    sol=solve_ivp(rhs,(zh,zz),[f0,P0],rtol=1e-11,atol=1e-14)
    f,P=sol.y[:,-1]; fpz0=P/muf(zz)
    out.append((fpz0>0, f>0, round(tt,3),round(x0/(tt-1),3),round(bb/(1+AA/zz),5)))
print(sum(o[0] for o in out),len(out)); print([o for o in out if not o[0]][:8])
