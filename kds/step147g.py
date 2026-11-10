import numpy as np, sympy as sp
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s_=z0-z; mus=z*(z-1)*(t-z); Ls=A/z+B-3*z; F=mus/s_**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mus*sp.diff(F,z),z)-2*s_*sp.diff(Ls*F,z))))
Nf=sp.lambdify((t,z0,A,B,S),sp.expand(nu.subs(z,z0-S)))
rng=np.random.default_rng(12); n=600000
x0=10**rng.uniform(-4,4,n); zz=1+x0
AA=10**rng.uniform(-5,8,n); gam=10**rng.uniform(-6,5,n)*rng.integers(0,2,n)
L1=AA*(1-1/(2*zz))+5*zz-3+gam
Tmax=np.minimum(2*x0*L1,2*x0**2*(AA/zz+3)); feas=Tmax>x0
tt=1+x0+(Tmax-x0)*rng.uniform(0,1,n)**0.2      # bias toward tmax
u=rng.uniform(1e-7,1-1e-7,n); q=1+x0*u; s=zz-q
BB=3*zz-AA/zz+AA/(2*zz)+2*zz+gam
L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q); mup=(q-1)*(tt-q)+q*(tt-q)-q*(q-1); H=2*mu+s*mup
N=Nf(tt,zz,AA,BB,s)
for th in (1.3,1.5,1.7,1.9):
    g=th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp
    viol=feas&(N>=0)&((g<=0)|(H-th*s*s*L<0))
    print(th,'violations',viol.sum(),'(g)',(feas&(N>=0)&(g<=0)).sum(),'(start)',(feas&(N>=0)&(H-th*s*s*L<0)).sum())
import mpmath as mp
mp.mp.dps=60
Nexpr=sp.expand(nu.subs(z,z0-S)); Nm=sp.lambdify((t,z0,A,B,S),Nexpr,'mpmath')
th=1.5; g=th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp
idx=np.where(feas&(N>=0)&(g<=0))[0][:300]
real=0
for i in idx:
    v=[mp.mpf(float(x)) for x in (tt[i],zz[i],AA[i],BB[i],s[i])]
    if Nm(*v)>=0: real+=1
print('exact N>=0 among flagged',real,len(idx))
print('x0 range of flagged', np.percentile(x0[idx],[0,50,100]),'A',np.percentile(AA[idx],[0,50,100]))
