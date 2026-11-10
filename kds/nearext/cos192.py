import numpy as np, time, warnings; warnings.filterwarnings('ignore')
exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
rng=np.random.default_rng(9); t0=time.time(); R=[]; viol=0
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,3.5); frac=rng.uniform(1e-4,0.99999); m=int(rng.integers(1,60))
    try: t,z0,A,B,q1=pdq(r0,1.0,r2,m,frac)
    except Exception: continue
    Lt=A/t+B-3*t
    if Lt>0: continue
    Q2=q1*q1; K=1.05*(-Lt)*t*(t-1)/(2*(t-z0)**2); r=K/Q2; a=z0-1; T=t-z0
    if a*(1+a)>T*(1+2*a): viol+=1
    b=T*np.linspace(0,1,4001); c=T-b; z=1+a+b
    al=(z*z*(r*(2+4*a+b-5*c)-3)/(1+c*r)).max(); al2=t*t*max(0,r*(2+4*a+T)-3)
    be=2*K*t*(2+r*T)
    R.append((A/(max(al,0)+be), A/(al2+be)))
R=np.array(R); print(len(R),'mu(z0)<0 viol',viol,'min margin (alpha sup)',R[:,0].min(),'(crude alpha)',R[:,1].min(), np.percentile(R[:,1],[1,5,50]))
