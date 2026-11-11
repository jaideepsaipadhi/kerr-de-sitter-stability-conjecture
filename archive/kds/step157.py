import sympy as sp, numpy as np
z,t,z0,A,B=sp.symbols('z t z0 A B')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
Nf=sp.lambdify((t,z0,A,B,z),sp.expand(nu)); H=2*mu+s*sp.diff(mu,z); Hf=sp.lambdify((t,z0,z),H)
rng=np.random.default_rng(2); n=2000000
x0=10**rng.uniform(-4,4,n); zz=1+x0; tt=zz+x0*10**rng.uniform(-6,7,n)
a0=np.maximum(0,zz*((tt-1)/(2*x0**2)-3))
AA=a0+10**rng.uniform(-6,9,n)*rng.integers(0,2,n)*(1+a0); gam=10**rng.uniform(-6,8,n)*rng.integers(0,2,n)
ell=AA/(2*zz)+2*tt+gam; BB=ell-AA/zz+3*zz
R3=2*x0*(AA+BB-3)>tt-1
q=1+x0*rng.uniform(1e-7,1-1e-7,n); s_=zz-q
N=Nf(tt,zz,AA,BB,q); Hv=Hf(tt,zz,q); L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
for th in (1.0,0.75,0.5,0.3):
    g=th*(2-th)*s_*s_*L*L+mu*(th*L-(2-th)*s_*Lp)
    print(th,'g<=0:',(R3&(g<=0)).sum(),' start fail:',(R3&(N>=0)&(Hv-th*s_*s_*L<0)).sum())
th=1.0; bad=R3&(N>=0)&(Hv-th*s_*s_*L<0)
p=2*x0*(AA+BB-3)/(tt-1)
for nm,v in (('x0',x0),('p',p),('A/a0',AA/np.maximum(a0,1e-300)),('(t-1)/x0',(tt-1)/x0),('u',(q-1)/x0),('gam/t',gam/tt)):
    print(nm,np.percentile(v[bad],[0,10,50,90,100]))
