import sympy as sp, numpy as np, mpmath as mp, sys
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
NP=sp.Poly(sp.expand(nu.subs(z,z0-S)),S); cf=[sp.lambdify((t,z0,A,B),c) for c in NP.all_coeffs()]
mp.mp.dps=30; CC=float(sys.argv[1]); DD=float(sys.argv[2]) if len(sys.argv)>2 else 0.0
rng=np.random.default_rng(4); tot=0; fail=0; worst=1e9
for i in range(3000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=CC*AA/zz+DD*zz+10**rng.uniform(-4,1); BB=3*zz-AA/zz+bb
    if 2*x0*(AA+BB-3)<=tt-1: continue
    c=[mp.mpf(g(tt,zz,AA,BB)) for g in cf]
    r=[mp.re(x) for x in mp.polyroots(c,maxsteps=300,extraprec=300) if abs(mp.im(x))<1e-15 and 0<mp.re(x)<x0]
    Sh=float(r[0]); zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1)
    H=2*mu_+Sh*mup; lam=H/Sh**2
    q=np.linspace(zh,zz,2001); sq=zz-q; Lq=AA/q+BB-3*q; Lp=-AA/q**2-3; muq=q*(q-1)*(tt-q)
    cond=lam**2*sq**2-2*lam*sq**2*Lq+muq*(2*sq*(-Lp)-lam)
    tot+=1; m=np.max(cond/(muq*lam+lam**2*sq**2+1e-300))
    worst=min(worst,-m)
    if m>0: fail+=1
print(CC,'tot',tot,'fail',fail,'best-worst margin',worst)
