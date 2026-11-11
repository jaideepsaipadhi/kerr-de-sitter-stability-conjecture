import numpy as np, time, sympy as sp, warnings; warnings.filterwarnings('ignore')
exec(open('phys149.py').read().split("rng=np.random.default_rng(1)")[0])
z_,t_,z0_,A_,B_=sp.symbols('z t z0 A B')
s=z0_-z_; mu=z_*(z_-1)*(t_-z_); L=A_/z_+B_-3*z_; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z_),z_)-2*s*sp.diff(L*F,z_))))
N=sp.expand(nu); H=2*mu+s*sp.diff(mu,z_); T=sp.expand(H-s**2*L)
NB=N.coeff(B_); Ts=sp.together(T+s**2*N/NB); n,d=sp.fraction(sp.factor(Ts)); n=sp.expand(n)
AT=sp.lambdify((t_,z0_,z_),-n.subs(A_,0)/n.coeff(A_))
Nf=sp.lambdify((t_,z0_,A_,B_,z_),N)
rng=np.random.default_rng(3); t0=time.time(); R=[]
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=10**rng.uniform(0.02,3.5); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,40))
    try: tt,zz,AA,BB=pdata(r0,1.0,r2,m,frac)
    except Exception: continue
    x0=zz-1
    if not (2*x0*(AA+BB-3)>tt-1 and 2*x0**2*(AA/zz+3)>=tt-1): continue
    q=1+x0*np.linspace(1e-5,1-1e-3,800); at=AT(tt,zz,q)
    i=np.where(at>=AA)[0]
    if len(i)==0: R.append((np.inf,x0,tt)); continue
    zT=q[i.min()]; Nv=Nf(tt,zz,AA,BB,zT)
    R.append((Nv,x0,tt))
R=np.array(R); print(len(R),'N(zT;Bc)>0 (violations):',(R[:,0]>0).sum(), ' zT exists frac',np.isfinite(R[:,0]).mean())
