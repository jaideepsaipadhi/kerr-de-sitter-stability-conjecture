import numpy as np, time, warnings; warnings.filterwarnings('ignore')
exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
rng=np.random.default_rng(5); t0=time.time(); R=[]
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,3.5); frac=rng.uniform(1e-4,0.99999); m=int(rng.integers(1,60))
    try: t,z0,A,B,q1=pdq(r0,1.0,r2,m,frac)
    except Exception: continue
    Lt=A/t+B-3*t
    if Lt>0: continue
    Q2=q1*q1; K=1.05*(-Lt)*t*(t-1)/(2*(t-z0)**2); r=K/Q2
    a=z0-1; T=t-z0; b=T*np.linspace(1e-6,1-1e-6,4000); c=T-b; z=1+a+b; z1=a+b
    g=z*(z-1)-c*(5*z-3*z0-1); p=2*a*a+2*a*b+2*a+b
    h2=a*a*b-6*a*a*c+2*a*b*b-6*a*b*c+a*b-6*a*c+b**3+b*b-3*b*c
    h3=a*a*b-2*a*a*c+2*a*b*b-2*a*b*c+a*b-2*a*c+b**3+b*b-b*c
    den=b*z1**2*(1+c*r)
    P1=(-3*b*z*z*z1**2+z*z*z1**2*g*r)/den; P2=b*b*r*(2*p-h2*r-c*h3*r*r)/den
    ex=(P1+Q2*P2).max(); cr=P1.max()+Q2*max(P2.max(),0)
    rc=0.7/T; denc=b*z1**2*(1+c*rc); P1c=(-3*b*z*z*z1**2+z*z*z1**2*g*rc)/denc; P2c=b*b*rc*(2*p-h2*rc-c*h3*rc*rc)/denc; cc=P1c.max()+Q2*max(P2c.max(),0)
    
    R.append((A/ex if ex>0 else np.inf, A/cr if cr>0 else np.inf, r*T, A/cc if cc>0 else np.inf))
R=np.array(R); f=np.isfinite(R[:,0]); g2=np.isfinite(R[:,1])
print(len(R),'exact min',R[f,0].min(),'crude min',R[g2,1].min(),'rT max',R[:,2].max())
h=np.isfinite(R[:,3]); print('capped min',R[h,3].min(), np.percentile(R[h,3],[1,5]))
