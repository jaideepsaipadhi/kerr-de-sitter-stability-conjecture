import numpy as np
rng=np.random.default_rng(8); worst=-9; arg=None
for _ in range(400000):
    t=1+10**rng.uniform(-2,5); z0=1+rng.uniform(0,1)*(t-1); z=1+rng.uniform(0,1)*(z0-1)
    B=5.5*t-1.25
    if z>=B/6: continue
    x=z-1; x0=z0-1; s=z0-z; mu0=z0*x0*(t-z0); d=B-6*z
    lhs=4.207*x*(t-z)*d/(B-3*z); rhs=np.sqrt(1.25*mu0*(B-3*z))
    v=(lhs-rhs)/(lhs+rhs)
    if v>worst: worst=v; arg=(t,(z0-1)/(t-1),(z-1)/x0)
print(worst,arg)
