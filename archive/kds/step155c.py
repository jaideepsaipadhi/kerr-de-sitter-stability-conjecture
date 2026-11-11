import sympy as sp, numpy as np
z,t,z0,A,B=sp.symbols('z t z0 A B')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
Nf=sp.lambdify((t,z0,A,B,z),sp.expand(nu)); H=2*mu+s*sp.diff(mu,z); Tf=sp.lambdify((t,z0,A,B,z),H-s**2*L)
rng=np.random.default_rng(2); n=2000000
x0=10**rng.uniform(-4,4,n); zz=1+x0; tt=zz+x0*10**rng.uniform(-6,7,n)
a0=np.maximum(0,zz*((tt-1)/(2*x0**2)-3))
AA=a0+10**rng.uniform(-6,9,n)*rng.integers(0,2,n)*(1+a0); gam=10**rng.uniform(-6,8,n)*rng.integers(0,2,n)
ell=AA/(2*zz)+2.5*tt+gam; BB=ell-AA/zz+3*zz
R3=2*x0*(AA+BB-3)>tt-1
q=1+x0*rng.uniform(1e-7,1-1e-7,n)
N=Nf(tt,zz,AA,BB,q); T=Tf(tt,zz,AA,BB,q)
bad=R3&(N>=0)&(T<0)
print('R3',R3.sum(),'bad',bad.sum())
if bad.any(): print('x0',np.percentile(x0[bad],[0,50,100]),'(t-z0)/x0',np.percentile(((tt-zz)/x0)[bad],[0,50,100]))
