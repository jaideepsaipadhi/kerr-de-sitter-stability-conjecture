import numpy as np
exec(open('step175.py').read().split("for cval")[0])
rng=np.random.default_rng(4); bad=[0,0]; n=0; ex=[]
while n<8000:
    x0=10**rng.uniform(-4,3); z0=1+x0; t=z0+x0*10**rng.uniform(-3,6)
    if (-3*z0*z0+2*(t+1)*z0-t)<0.5*z0*t: continue
    q1=10**rng.uniform(-3,4); A=1e-3*q1*q1*z0*z0/t*10**rng.uniform(0,6)
    kmax=min(0.6*q1/(t-z0), 0.2*q1*(t-z0)*(A/t**2+3)/(t*(t-1)))
    kap=kmax*rng.uniform(0,1)
    a,b=check(t,z0,A,q1,kap); bad[0]+=not a; bad[1]+=not b; n+=1
    if not b and len(ex)<5: ex.append((x0,t/z0,A,q1,kap))
print('Cfail',bad[0],'Dfail',bad[1],'of',n); print(ex)
