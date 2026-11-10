import numpy as np
exec(open('step175.py').read().split("for cval")[0])
rng=np.random.default_rng(3); bad=[0,0]; n=0
while n<8000:
    x0=10**rng.uniform(-4,3); z0=1+x0; t=z0+x0*10**rng.uniform(-3,4)
    A=10**rng.uniform(-3,8)*rng.integers(0,2); q1=10**rng.uniform(-3,4)
    kmax=min(0.6*q1/(t-z0), 0.2*q1*(t-z0)*(A/t**2+3)/(t*(t-1)))
    kap=kmax*rng.uniform(0,1)
    a,b=check(t,z0,A,q1,kap); bad[0]+=not a; bad[1]+=not b; n+=1
print('Cfail',bad[0],'Dfail',bad[1],'of',n)
rng=np.random.default_rng(3); k=0
while k<6:
    x0=10**rng.uniform(-4,3); z0=1+x0; t=z0+x0*10**rng.uniform(-3,4)
    A=10**rng.uniform(-3,8)*rng.integers(0,2); q1=10**rng.uniform(-3,4)
    kmax=min(0.6*q1/(t-z0), 0.2*q1*(t-z0)*(A/t**2+3)/(t*(t-1)))
    kap=kmax*rng.uniform(0,1)
    z=z0+(t-z0)*np.linspace(1e-6,1-1e-6,600); mu=z*(z-1)*(t-z); mup=(z-1)*(t-z)+z*(t-z)-z*(z-1)
    Q=q1+kap*(t-z); qt=(z-z0)*Q
    c=lambda zz: kap*(zz-z0)**2*(2*q1+kap*(t-zz))/(zz*(zz-1)); h=1e-7*(t-z0); dc=(c(z+h)-c(z-h))/(2*h)
    D=qt*(-A/z**2-3+dc)-kap*mup
    if D.max()>0:
        i=D.argmax(); k+=1
        print('u=%.3f x0=%.3g t/z0=%.3g A=%.3g q1=%.3g | qtL\'=%.3g qt*dc=%.3g -kap mu\'=%.3g'%((z[i]-z0)/(t-z0),x0,t/z0,A,q1,qt[i]*(-A/z[i]**2-3),qt[i]*dc[i],-kap*mup[i]))
