import numpy as np, sympy as sp, sys
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s_=z0-z; mus=z*(z-1)*(t-z); Ls=A/z+B-3*z; F=mus/s_**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mus*sp.diff(F,z),z)-2*s_*sp.diff(Ls*F,z))))
Nf=sp.lambdify((t,z0,A,B,S),sp.expand(nu.subs(z,z0-S)))
rng=np.random.default_rng(12); n=600000
x0=10**rng.uniform(-4,4,n); zz=1+x0
AA=10**rng.uniform(-5,8,n); gam=10**rng.uniform(-6,5,n)*rng.integers(0,2,n)
BB=3*zz-AA/zz+AA/(2*zz)+2*zz+gam
L1=AA+BB-3
Tmax=np.minimum(2*x0*L1,2*x0**2*(AA/zz+3)); feas=Tmax>x0
tt=1+x0+(Tmax-x0)*rng.uniform(0,1,n)**0.2
u=rng.uniform(1e-7,1-1e-7,n); q=1+x0*u; s=zz-q
L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q); mup=(q-1)*(tt-q)+q*(tt-q)-q*(q-1); H=2*mu+s*mup
N=Nf(tt,zz,AA,BB,s)
best=np.zeros(n,bool)
for th in np.linspace(0.3,1.95,34):
    g=th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp
    best|= (g>0)&(H-th*s*s*L>=0)
bad=feas&(N>=0)&~best
print('bad',bad.sum())
for nm,c in (('B>=0',BB>=0),('B>=t',BB>=tt),('B>=-A/(4z0)',BB>=-AA/(4*zz)),('beta>=0.6A/z0',(AA/zz+BB-3*zz)>=0.6*AA/zz)):
    print(nm,'bad remaining',(bad&c).sum(),'feasible pts',(feas&(N>=0)&c).sum())
r=BB/(AA/zz)
print('bad B/(A/z0)',np.percentile(r[bad],[0,50,100]),' bad (B-t)/t',np.percentile(((BB-tt)/tt)[bad],[0,50,100]),' bad x0',np.percentile(x0[bad],[0,50,100]),'bad t/x0',np.percentile((tt/x0)[bad],[0,50,100]))
