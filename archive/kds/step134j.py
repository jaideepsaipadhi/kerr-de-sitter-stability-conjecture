import sympy as sp, numpy as np, mpmath as mp
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
N=sp.Poly(sp.expand(nu.subs(z,z0-S)),S)
mp.mp.dps=50
for tt,zz in ((1.019,1.007),(2.102,1.06),(1.018,1.017)):
  tt=mp.mpf(tt); zz=mp.mpf(zz); x0=zz-1
  a0=zz*((tt-1)/(2*x0**2)-3)
  for fa,bb in ((1,0),(1.01,0),(1.01,0.01),(1.1,0.1),(1,1)):
    AA=a0*fa; BB=3*zz-AA/zz+bb
    p=2*x0*(AA+BB-3)/(tt-1)
    cs=[mp.mpf(str(sp.N(c.subs({t:tt,z0:zz,A:AA,B:BB}),60))) for c in N.all_coeffs()]
    rts=[r for r in mp.polyroots(cs,maxsteps=200,extraprec=200) if abs(mp.im(r))<1e-20 and 0<mp.re(r)<x0]
    print(float(tt),float(zz),fa,bb,'p=%.4f'%float(p),'roots',len(rts),'N(x0) sign',mp.sign(mp.polyval(cs,x0)))
