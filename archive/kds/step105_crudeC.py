import numpy as np
rng=np.random.default_rng(0); worst=-9; arg=None
for _ in range(400000):
    t=1+10**rng.uniform(-2,5); zf=rng.uniform(0,1); z0=1+zf*(t-1)
    B=(5.5*t-1.25)*10**rng.uniform(0,3); A=B*10**rng.uniform(-2,6)
    y=rng.uniform(0,1); z=1+y*(z0-1); s=z0-z
    L=A/z+B-3*z
    if L<=0: continue
    if 2*(z0-1)*(A+B-3)/(t-1)<=1: continue   # R3 only
    mu=z*(z-1)*(t-z); mup=-3*z*z+2*(1+t)*z-t; H=2*mu+s*mup
    g=(A/z**2+3)/L
    Cc=s*mu*(1.7*g+0.45/((z-1)*z))-H
    m=Cc/(abs(s*mu*(1.7*g+0.45/((z-1)*z)))+abs(H))
    if m>worst: worst=m; arg=(t,zf,A/B,y)
print(worst,arg)
