import sympy as sp, numpy as np, mpmath as mp
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
NP=sp.Poly(sp.expand(nu.subs(z,z0-S)),S); cf=[sp.lambdify((t,z0,A,B),c) for c in NP.all_coeffs()]
mp.mp.dps=30
rng=np.random.default_rng(4); rows=[]
for i in range(3000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=10**rng.uniform(-4,2)*(1+AA/zz); BB=3*zz-AA/zz+bb
    p=2*x0*(AA+BB-3)/(tt-1)
    if p<=1: continue
    c=[mp.mpf(g(tt,zz,AA,BB)) for g in cf]
    r=[mp.re(x) for x in mp.polyroots(c,maxsteps=300,extraprec=300) if abs(mp.im(x))<1e-15 and 0<mp.re(x)<x0]
    Sh=float(r[0]); zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1)
    H=2*mu_+Sh*mup; Q=H/Sh-2*(AA*(zz/zh-1-np.log(zz/zh))+1.5*Sh**2)
    rows.append((Q>0,p,bb/(1+AA/zz),AA/(a0+1e-12) if a0>0 else np.inf,Sh/x0,tt,x0/(tt-1)))
R=np.array([r[1:] for r in rows if not r[0]]); G=np.array([r[1:] for r in rows if r[0]])
for nm,i in (('p',0),('beta_n',1),('A/a0',2),('Sh/x0',3),('t',4),('zf',5)):
    print(nm,'fail range',np.round(np.percentile(R[:,i],[0,50,100]),4),' pass range',np.round(np.percentile(G[:,i],[0,50,100]),4))
