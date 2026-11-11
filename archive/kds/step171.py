import numpy as np, sympy as sp, warnings; warnings.filterwarnings('ignore')
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s_=z0-z; mus=z*(z-1)*(t-z); Ls=A/z+B-3*z; F=mus/s_**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mus*sp.diff(F,z),z)-2*s_*sp.diff(Ls*F,z))))
Nf=sp.lambdify((t,z0,A,B,S),sp.expand(nu.subs(z,z0-S)))
rng=np.random.default_rng(16); n=400000
x0=10**rng.uniform(-4,4,n); zz=1+x0; tt=zz+x0*10**rng.uniform(-6,7,n)
AA=10**rng.uniform(-6,9,n)*rng.integers(0,2,n)
Bm=AA/(2*zz)+2*tt-AA/zz+3*zz
Sh=x0*rng.uniform(1e-6,1-1e-6,n)
gate=Nf(tt,zz,AA,Bm,Sh)>=0
zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); w0=(2*mu_+Sh*mup)/Sh
q=zh+Sh*rng.uniform(0,1-1e-9,n); s=zz-q; L=AA/q+Bm-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
Q=w0*L-2*mu*Lp
print('gate pts',gate.sum(),'Q<=0 among gate',(gate&(Q<=0)).sum(),' Q<=0 overall',(Q<=0).sum())
