import numpy as np, time, warnings; warnings.filterwarnings('ignore')
exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
rng=np.random.default_rng(9); t0=time.time(); R=[]
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,3.5); frac=rng.uniform(1e-4,0.99999); m=int(rng.integers(1,60))
    try: t,z0,A,B,q1=pdq(r0,1.0,r2,m,frac)
    except Exception: continue
    Lt=A/t+B-3*t
    if Lt>0: continue
    Q2=q1*q1; K=1.05*(-Lt)*t*(t-1)/(2*(t-z0)**2); r=K/Q2; a=z0-1; T=t-z0
    b=T*np.linspace(1e-6,1,4001); c=T-b; z=1+a+b; z1=a+b
    g=z*(z-1)-c*(1+2*a+5*b); p=2*a*a+2*a*b+2*a+b
    h2=a*a*b-6*a*a*c+2*a*b*b-6*a*b*c+a*b-6*a*c+b**3+b*b-3*b*c
    h3=a*a*b-2*a*a*c+2*a*b*b-2*a*b*c+a*b-2*a*c+b**3+b*b-b*c
    P1=z*z*(g*r-3*b)/(b*(1+c*r)); P2=b*r*(2*p-h2*r-c*h3*r*r)/(z1**2*(1+c*r))
    al=(z*z*(r*(2+4*a+b-5*c)-3)/(1+c*r)).max()
    R.append((max(P1.max(),1e-300), max(al,1e-300), max(Q2*P2.max(),1e-300), 2*K*t*(2+r*T), A))
R=np.array(R); print('alpha bound/true  max',np.max(R[:,1]/R[:,0]),' A/alpha_bd min',np.min(R[:,4]/R[:,1]))
print('beta bound/true max',np.max(R[:,3]/R[:,2]),' A/beta_bd min',np.min(R[:,4]/R[:,3]), 'A/beta_true min',np.min(R[:,4]/R[:,2]))
