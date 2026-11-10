import sympy as sp, numpy as np, mpmath as mp
from scipy.integrate import solve_ivp
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
NP=sp.Poly(sp.expand(nu.subs(z,z0-S)),S); cf=[sp.lambdify((t,z0,A,B),c) for c in NP.all_coeffs()]
mp.mp.dps=30
import sys; CC=float(sys.argv[1])
rng=np.random.default_rng(4); res=[]
for i in range(3000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=CC*AA/zz+10**rng.uniform(-4,2)*(1+AA/zz)*0.0+10**rng.uniform(-4,1); BB=3*zz-AA/zz+bb
    if 2*x0*(AA+BB-3)<=tt-1: continue
    c=[mp.mpf(g(tt,zz,AA,BB)) for g in cf]
    r=[mp.re(x) for x in mp.polyroots(c,maxsteps=300,extraprec=300) if abs(mp.im(x))<1e-15 and 0<mp.re(x)<x0]
    Sh=float(r[0]); zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1)
    H=2*mu_+Sh*mup; Q=H/Sh-2*(AA*(zz/zh-1-np.log(zz/zh))+1.5*Sh**2)
    if Q>0: res.append((True,True,9)); continue
    print("Qfail",round(tt,3),round(x0/(tt-1),3),AA,bb)
    Lf=lambda q: AA/q+BB-3*q; muf=lambda q: q*(q-1)*(tt-q)
    f0=1.0; Phi0=f0*(H/Sh-2*Sh*Lf(zh))
    def rhs(q,y):  # frozen-Phi lower bound f_l and accumulated integral
        fl,I=y; return [(Phi0+2*(zz-q)*Lf(q)*fl)/muf(q), 2*Lf(q)*fl]
    sol=solve_ivp(rhs,(zh,zz),[f0,0.0],rtol=1e-10,atol=1e-14)
    fl,I=sol.y[:,-1]; minfl=sol.y[0].min()
    res.append((Phi0+I>0, minfl>0, (Phi0+I)/abs(Phi0)))
print(len(res),'pass',sum(r[0] and r[1] for r in res),'min ratio',min(r[2] for r in res))
