import sympy as sp, numpy as np, warnings; warnings.filterwarnings('ignore')
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s_=z0-z; mus=z*(z-1)*(t-z); Ls=A/z+B-3*z; F=mus/s_**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mus*sp.diff(F,z),z)-2*s_*sp.diff(Ls*F,z))))
NP=sp.Poly(sp.expand(nu.subs(z,z0-S)),S); cf=[sp.lambdify((t,z0,A,B),c) for c in NP.all_coeffs()]
import mpmath as mp
rng=np.random.default_rng(3); tot=startfail=altfail=0
import time; T0=time.time()
while time.time()-T0<45:
    x0=10**rng.uniform(-4,4); zz=1+x0; tt=zz+x0*10**rng.uniform(-6,7)
    a0=max(0,zz*((tt-1)/(2*x0**2)-3)); AA=a0+10**rng.uniform(-6,9)*rng.integers(0,2)*(1+a0); gam=10**rng.uniform(-6,8)*rng.integers(0,2)
    ell=AA/(2*zz)+2*tt+gam; BB=ell-AA/zz+3*zz
    pass
    c=[mp.mpf(g(tt,zz,AA,BB)) for g in cf]
    try: r=[mp.re(v) for v in mp.polyroots(c,maxsteps=300,extraprec=100) if abs(mp.im(v))<1e-25 and 0<mp.re(v)<x0]
    except Exception: continue
    if len(r)!=1: continue
    tot+=1; Sh=float(r[0]); zh=zz-Sh
    mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); w0=(2*mu_+Sh*mup)/Sh; Lh=AA/zh+BB-3*zh
    if w0>=Sh*Lh: continue
    startfail+=1
    q=np.linspace(zh,zz,4000); s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    seg=s*L>=w0
    R=-w0*w0+2*s*L*w0-2*mu*s*Lp
    altfail+= np.any(seg&(R<=0))
print('tot',tot,'start fail',startfail,'alt (constant-w0 then sL) fail',altfail)
