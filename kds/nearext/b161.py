import numpy as np, time, sympy as sp
exec(open('phys149.py').read().split("rng=np.random.default_rng(1)")[0])
zz_,tt_,z0_,A_,B_=sp.symbols('z t z0 A B')
s=z0_-zz_; mu=zz_*(zz_-1)*(tt_-zz_); L=A_/zz_+B_-3*zz_; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,zz_),zz_)-2*s*sp.diff(L*F,zz_))))
N=sp.expand(nu); H=2*mu+s*sp.diff(mu,zz_); T=sp.expand(H-s**2*L)
Bstar=sp.solve(N,B_)[0]; Tstar=sp.together(T.subs(B_,Bstar))
f=sp.lambdify((tt_,z0_,A_,zz_),[Tstar,Bstar])
FAILS=[]; UP=[]; LO=[]; import warnings; warnings.filterwarnings('ignore'); rng=np.random.default_rng(12); t0=time.time(); cs=(1.0,1.5,2,4,10); ok={c:0 for c in cs}; n=0
while time.time()-t0<45:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=10**rng.uniform(0.02,3.5); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,40))
    try: tt,zz,AA,BB=pdata(r0,1.0,r2,m,frac)
    except Exception: continue
    x0=zz-1
    if not (2*x0*(AA+BB-3)>tt-1 and 2*x0**2*(AA/zz+3)>=tt-1): continue
    q=1+x0*np.linspace(1e-6,1-1e-4,600); Ts,Bs=f(tt,zz,AA,q); n+=1
    zs=1+0.62*x0*x0/zz; sel=(q>=zs)&np.isfinite(Ts)
    UP.append(Ts[sel].min()/np.abs(Ts[sel]).max() if sel.any() else 1); mono=np.all(np.diff(Bs[np.isfinite(Bs)])>=-1e-9*np.abs(Bs[np.isfinite(Bs)][1:]))
    LO.append(0 if mono else 1)
print(n,ok)

for F_ in FAILS[:8]: print(F_)
import pickle; pickle.dump(FAILS,open('fails159.pkl','wb'))

UP=np.array(UP); LO=np.array(LO); print('pts',len(UP),'T* min normalised on [z*,z0)',UP.min(),'#neg',(UP<0).sum(),' B* nonmonotone pts',LO.sum())
