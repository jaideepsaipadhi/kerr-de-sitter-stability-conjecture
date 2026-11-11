import numpy as np, time, warnings; warnings.filterwarnings('ignore')
exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
rng=np.random.default_rng(11); t0=time.time(); P=[0,1,2,4]; bad={p:[0,0] for p in P}; marg={p:[] for p in P}; n=0
while time.time()-t0<45:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,3.5); frac=rng.uniform(1e-4,0.99999); m=int(rng.integers(1,60))
    try: t,z0,A,B,q1=pdq(r0,1.0,r2,m,frac)
    except Exception: continue
    Lt=A/t+B-3*t
    if Lt>0: continue
    n+=1; kap=1.05*(-Lt)*t*(t-1)/(2*q1*(t-z0)**2)
    z=z0+(t-z0)*np.linspace(1e-5,1-1e-5,4000); mu=z*(z-1)*(t-z)
    L=A/z+B-3*z; q=q1*(z-z0)
    for p in P:
        phi=((z-z0)/(t-z0))**p; qt=(z-z0)*(q1+kap*(t-z)*phi)
        Lt_=L+(qt**2-q**2)/mu
        d=lambda f: np.gradient(f,z)
        C=d(mu)*qt-2*mu*d(qt); D=qt*d(Lt_)+0.5*d(mu*d(d(qt)))
        s=slice(20,-20)
        bad[p][0]+= (C[s]>1e-9*abs(C[s]).max()).any(); bad[p][1]+= (D[s]>1e-9*abs(D[s]).max()).any()
print('gate pts',n, {p:bad[p] for p in P})
