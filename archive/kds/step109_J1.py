import numpy as np
rng=np.random.default_rng(7); W={}; 
for key in ('A0+bulk','A0+bulk+R3'): W[key]=(-9,None)
for _ in range(400000):
    t=1+10**rng.uniform(-2,5); z0=1+rng.uniform(0,1)*(t-1); z=1+rng.uniform(0,1)*(z0-1)
    B=5.5*t-1.25
    if z>=B/6: continue
    x=z-1; x0=z0-1; s=z0-z; mu0=z0*x0*(t-z0); d=B-6*z
    Ab=z0*(10*mu0/(8*s*s)-B+3*z0)
    Ar=max((t-1)/(2*x0)-B+3, z0*((t-1)/(2*x0*x0)-3))
    for key,Am in (('A0+bulk',max(0,Ab)),('A0+bulk+R3',max(0,Ab,Ar))):
        L=Am/z+B-3*z
        lhs=4.207*x*(t-z)*d/s; v=(lhs-L*L)/(lhs+L*L)
        if v>W[key][0]: W[key]=(v,(t,(z0-1)/(t-1),(z-1)/x0))
print(W)
