import numpy as np
rng=np.random.default_rng(2)
def check(t,z0,A,q1,kap,n=600):
    z=z0+(t-z0)*np.linspace(1e-6,1-1e-6,n); mu=z*(z-1)*(t-z); mup=(z-1)*(t-z)+z*(t-z)-z*(z-1)
    Q=q1+kap*(t-z); qt=(z-z0)*Q; qtp=Q-kap*(z-z0)
    c=lambda zz: kap*(zz-z0)**2*(2*q1+kap*(t-zz))/(zz*(zz-1))
    h=1e-7*(t-z0); dc=(c(z+h)-c(z-h))/(2*h)
    Ltp=-A/z**2-3+dc
    C=mup*qt-2*mu*qtp; D=qt*Ltp-kap*mup
    return C.max()<=0, D.max()<=0
for cval in (0.1,0.3,1.0,3.0):
    bad=[0,0]; n=0
    for _ in range(6000):
        x0=10**rng.uniform(-3,3); z0=1+x0; t=z0+x0*10**rng.uniform(-3,3)
        A=10**rng.uniform(-3,6)*rng.integers(0,2); q1=10**rng.uniform(-3,4)
        kap=cval*q1/(t-z0)*rng.uniform(0,1)
        a,b=check(t,z0,A,q1,kap); bad[0]+=not a; bad[1]+=not b; n+=1
    print('kappa<= c q1/(t-z0), c=',cval,'Cfail',bad[0],'Dfail',bad[1],'of',n)
