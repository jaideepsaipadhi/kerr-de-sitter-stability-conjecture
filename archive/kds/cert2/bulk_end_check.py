import numpy as np
rng=np.random.default_rng(21); cnt={'convex':0,'notconvex':0,'endmax':0,'intmax':0,'n':0}
def B3(t,z0,z,A,B):
    x=z-1;x0=z0-1;s=z0-z;mu=z*x*(t-z);mup=-3*z*z+2*(1+t)*z-t
    L=A/z+B-3*z; g=(A/z**2+3)/L; v=1-x/x0; q=0.45*v/(x*z); qp=0.45*(-(1/x0)/(x*z)-v*(2*z-1)/(x*z)**2)
    h=np.where(mup>0,1.7*g+q,17/11*g+q)
    return 0.5*mup*h+0.85*mu*g*g+0.5*mu*(1.7*g+q)**2+0.5*mu*qp-0.5455*s*L*g-s*L*q
for _ in range(20000):
    t=1+10**rng.uniform(-2,5); z0=1+rng.uniform(0,1)*(t-1); z=1+rng.uniform(0,1)*(z0-1)
    x0=z0-1; s=z0-z; mu0=z0*x0*(t-z0); fl=5.5*t-1.25; L0s=1.25*mu0/s**2
    Bhi=L0s+3*z0
    if Bhi<=fl: continue
    Bs=np.linspace(fl,Bhi,400); As=z0*(L0s-Bs+3*z0)
    ok=(2*x0*(As+Bs-3)/(t-1)>1)&(2*x0**2*(As/z0+3)>=t-1)
    if ok.sum()<5: continue
    f=B3(t,z0,z,As,Bs)[ok]; cnt['n']+=1
    d2=np.diff(f,2)
    cnt['convex' if np.all(d2>=-1e-9*np.abs(f).max()) else 'notconvex']+=1
    i=np.argmax(f); cnt['endmax' if i in (0,len(f)-1) else 'intmax']+=1
print(cnt)
