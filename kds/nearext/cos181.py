import numpy as np, time, warnings; warnings.filterwarnings('ignore')
exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
rng=np.random.default_rng(7); t0=time.time(); R=[]
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,3.5); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,40))
    try: tt,zz,AA,BB,q1=pdq(r0,1.0,r2,m,frac)
    except Exception: continue
    Lt=AA/tt+BB-3*tt
    if Lt>0: continue
    kap=1.05*(-Lt)*tt*(tt-1)/(2*q1*(tt-zz)**2)
    z=zz+(tt-zz)*np.linspace(1e-6,1-1e-6,3000); mu=z*(z-1)*(tt-z); mup=(z-1)*(tt-z)+z*(tt-z)-z*(z-1)
    Q=q1+kap*(tt-z); qt=(z-zz)*Q
    cf=lambda w: kap*(w-zz)**2*(2*q1+kap*(tt-w))/(w*(w-1)); h=1e-7*(tt-zz); dc=(cf(z+h)-cf(z-h))/(2*h)
    num=qt*(dc-3)-kap*mup; coef=qt/z**2   # D = num - A*coef
    Areq=(num/coef).max()
    R.append((AA/Areq if Areq>0 else np.inf, AA*tt/(q1*zz)**2, kap/q1*(tt-zz), (zz-1)/(tt-1), tt))
R=np.array(R); fin=np.isfinite(R[:,0]); print(len(R),'A/Areq min',R[fin,0].min(),'pct',np.percentile(R[fin,0],[0,5,50]),'Areq<=0 frac',(~fin).mean())
